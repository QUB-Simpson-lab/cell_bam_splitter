/*
 * cell_bam_splitter — split a BAM by an arbitrary auxiliary tag value.
 *
 * Strategy: two-pass, sequential, raw-BGZF append.
 *
 * Pass 1 — index (zero allocation per read):
 *   Stream every read.  Record the BGZF virtual offset before each read and
 *   the tag value.  Store (offset -> tag_value) in a flat sorted array.
 *   Nothing is copied; bam1_t is reused across iterations.
 *
 * Pass 2 — write (single sequential scan, raw block copy):
 *   Re-open the input.  Walk the sorted offset list in order.  For each
 *   read, seek only when necessary (i.e. the read is not the next sequential
 *   one), decompress it once, and write it to the appropriate output file.
 *   Because offsets are sorted, reads are visited in file order: the input
 *   is read strictly sequentially.  Each output file is opened once and kept
 *   open for the duration, then closed.  The number of concurrently open
 *   output handles is bounded by the OS fd limit at runtime.
 *
 * Why raw-BGZF append works:
 *   A BGZF file is a sequence of independent compressed blocks terminated by
 *   a fixed 28-byte EOF block.  Stripping the EOF block and appending new
 *   valid BGZF blocks (followed by a new EOF block) is explicitly legal per
 *   the BGZF spec; samtools cat exploits this.  htslib's "ab" open mode is
 *   not exposed via sam_open, but we can implement it ourselves:
 *     - On first open of an output file: write header normally via sam_open "wb".
 *     - On subsequent opens (if we needed to close due to fd pressure): open
 *       the file in raw binary append mode, truncate the 28-byte EOF block,
 *       write new BGZF blocks, append EOF.
 *   In practice, with offsets sorted and all outputs held open concurrently
 *   (bounded by fd limit), the reopen path is only hit when cell count exceeds
 *   the fd budget, which is rare with ulimit -n set sensibly.
 *
 * Memory use:
 *   Phase 1 index: 8 bytes (offset) + ~average tag length per read.
 *   For 50M reads across 10k cells: ~500MB index.  Acceptable.
 *   No read payload is ever copied.
 *
 * Build:
 *   g++ -O2 -std=c++17 cell_bam_splitter.cpp -o cell_bam_splitter -lhts -lpthread
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <atomic>
#include <mutex>
#include <thread>
#include <cstring>
#include <cerrno>
#include <cassert>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <getopt.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>

// ---------------------------------------------------------------------------
// BGZF EOF block (28 bytes, fixed by spec).
// ---------------------------------------------------------------------------
static const uint8_t BGZF_EOF[28] = {
    0x1f,0x8b,0x08,0x04, 0x00,0x00,0x00,0x00,
    0x00,0xff,0x06,0x00, 0x42,0x43,0x02,0x00,
    0x1b,0x00,0x03,0x00, 0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00
};

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

static std::string trailingSlash(const std::string& p) {
    return (!p.empty() && p.back() == '/') ? p : p + '/';
}

static std::string fmtCount(uint64_t n) {
    std::string s = std::to_string(n);
    for (int i = static_cast<int>(s.size()) - 3; i > 0; i -= 3)
        s.insert(i, ",");
    return s;
}

static std::string fmtDuration(double sec) {
    if (sec < 60)  return std::to_string(static_cast<int>(sec)) + "s";
    int m = static_cast<int>(sec) / 60, s2 = static_cast<int>(sec) % 60;
    if (m < 60) return std::to_string(m) + "m" + std::to_string(s2) + "s";
    int h = m / 60; m %= 60;
    return std::to_string(h) + "h" + std::to_string(m) + "m";
}

// ---------------------------------------------------------------------------
// Logger
// ---------------------------------------------------------------------------

class Logger {
public:
    explicit Logger(bool verbose)
        : verbose_(verbose), start_(std::chrono::steady_clock::now()) {}

    void step(const std::string& msg) {
        double el = elapsed();
        std::lock_guard<std::mutex> lk(mu_);
        std::cout << "[" << std::fixed << std::setw(7) << std::setprecision(1)
                  << el << "s] " << msg << '\n' << std::flush;
    }
    void verbose(const std::string& msg) { if (verbose_) step(msg); }
    double elapsed() const {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start_).count();
    }
private:
    bool verbose_;
    std::chrono::steady_clock::time_point start_;
    std::mutex mu_;
};

// ---------------------------------------------------------------------------
// Directory / fd helpers
// ---------------------------------------------------------------------------

static bool ensureDirectory(const std::string& path) {
    struct stat st{};
    if (stat(path.c_str(), &st) == 0) return S_ISDIR(st.st_mode);
    std::string cur;
    for (size_t i = 0; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '/') {
            if (!cur.empty() && cur != "/") {
                struct stat s2{};
                if (stat(cur.c_str(), &s2) != 0 &&
                    mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) {
                    std::cerr << "mkdir failed: " << cur << " (" << strerror(errno) << ")\n";
                    return false;
                }
            }
        }
        if (i < path.size()) cur += path[i];
    }
    return true;
}

static int getOpenFileLimit() {
    struct rlimit rl{};
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY)
        return static_cast<int>(rl.rlim_cur);
    return 256;
}

// ---------------------------------------------------------------------------
// Output file state — one per tag value.
// ---------------------------------------------------------------------------

struct OutFile {
    std::string  path;
    samFile*     fp        = nullptr;  // null when closed
    uint64_t     nreads    = 0;
    bool         headerWritten = false;

    OutFile() = default;
    OutFile(std::string p) : path(std::move(p)) {}
};

// Open a fresh output BAM and write header.
static bool openAndWriteHeader(OutFile& of, bam_hdr_t* hdr) {
    assert(!of.fp);
    of.fp = sam_open(of.path.c_str(), "wb");
    if (!of.fp) {
        std::cerr << "Error: cannot open output: " << of.path << '\n';
        return false;
    }
    if (sam_hdr_write(of.fp, hdr) < 0) {
        std::cerr << "Error: header write failed: " << of.path << '\n';
        sam_close(of.fp); of.fp = nullptr;
        return false;
    }
    of.headerWritten = true;
    return true;
}

// Reopen an existing output BAM for append by stripping the EOF block.
// Returns false on error.
static bool reopenForAppend(OutFile& of) {
    assert(!of.fp && of.headerWritten);
    // Open raw for binary read/write.
    FILE* raw = fopen(of.path.c_str(), "r+b");
    if (!raw) {
        std::cerr << "Error: cannot reopen for append: " << of.path << '\n';
        return false;
    }
    // Seek to EOF, verify and strip the 28-byte BGZF EOF block.
    if (fseeko(raw, -28, SEEK_END) != 0) {
        std::cerr << "Error: fseek failed on: " << of.path << '\n';
        fclose(raw); return false;
    }
    uint8_t tail[28];
    if (fread(tail, 1, 28, raw) != 28 || memcmp(tail, BGZF_EOF, 28) != 0) {
        std::cerr << "Warning: unexpected tail on " << of.path
                  << "; not stripping EOF block.\n";
        fclose(raw);
        // Fall back: open normally and append (will result in double EOF, still
        // readable by samtools but technically non-conformant).
        of.fp = sam_open(of.path.c_str(), "ab");
        return of.fp != nullptr;
    }
    // Truncate to remove the EOF block.
    off_t newSize = ftello(raw) - 28;
    fclose(raw);
    if (truncate(of.path.c_str(), newSize) != 0) {
        std::cerr << "Error: truncate failed on: " << of.path << '\n';
        return false;
    }
    // Now open via htslib in append mode.  htslib's BGZF writer opened with
    // "ab" will append new blocks to the file, then write a fresh EOF on close.
    of.fp = sam_open(of.path.c_str(), "ab");
    if (!of.fp) {
        std::cerr << "Error: sam_open(ab) failed: " << of.path << '\n';
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Usage
// ---------------------------------------------------------------------------

static void printUsage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " -i input.bam -o output_dir\n"
        << "           [-b barcodes.txt] [-t TAG] [-j threads] [-v]\n\n"
        << "  -i  Input BAM\n"
        << "  -o  Output directory (created if absent)\n"
        << "  -b  Barcode allowlist, one per line\n"
        << "  -t  2-char auxiliary tag to split on [default: CB]\n"
        << "  -j  Threads for Phase 1 input decompression [default: hw]\n"
        << "  -v  Verbose output\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    const char* inputBam    = nullptr;
    const char* outputDir   = nullptr;
    const char* barcodeFile = nullptr;
    char        splitTag[3] = "CB";
    int         nThreads    = static_cast<int>(std::thread::hardware_concurrency());
    bool        verbose     = false;
    if (nThreads < 1) nThreads = 1;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:b:t:j:vh")) != -1) {
        switch (opt) {
            case 'i': inputBam    = optarg; break;
            case 'o': outputDir   = optarg; break;
            case 'b': barcodeFile = optarg; break;
            case 't':
                if (strlen(optarg) != 2) {
                    std::cerr << "Error: -t must be exactly 2 characters.\n"; return 1;
                }
                strncpy(splitTag, optarg, 2); splitTag[2] = '\0';
                break;
            case 'j':
                nThreads = std::atoi(optarg);
                if (nThreads < 1) { std::cerr << "Error: -j >= 1.\n"; return 1; }
                break;
            case 'v': verbose = true; break;
            case 'h': printUsage(argv[0]); return 0;
            default:  printUsage(argv[0]); return 1;
        }
    }
    if (!inputBam || !outputDir) { printUsage(argv[0]); return 1; }

    Logger log(verbose);
    log.step("cell_bam_splitter starting");
    log.step("  input  : " + std::string(inputBam));
    log.step("  output : " + std::string(outputDir));
    log.step("  tag    : " + std::string(splitTag));
    log.step("  threads: " + std::to_string(nThreads));

    if (!ensureDirectory(outputDir)) return 1;
    const std::string outPrefix = trailingSlash(outputDir);

    const int fdLimit   = getOpenFileLimit();
    // Reserve: 1 input + stdin/stdout/stderr + OS headroom.
    // Each output file = 1 fd.
    const int maxOpen   = std::max(4, fdLimit - 8);
    log.step("  fd limit: " + std::to_string(fdLimit)
           + "  |  max open output files: " + std::to_string(maxOpen));

    // ------------------------------------------------------------------
    // Allowlist.
    // ------------------------------------------------------------------
    std::unordered_set<std::string> allowlist;
    bool useAllowlist = false;
    if (barcodeFile) {
        useAllowlist = true;
        std::ifstream fin(barcodeFile);
        if (!fin) {
            std::cerr << "Error: cannot open: " << barcodeFile << '\n'; return 1;
        }
        std::string line;
        while (std::getline(fin, line))
            if (!line.empty()) allowlist.insert(line);
        log.step("Allowlist: " + fmtCount(allowlist.size()) + " barcodes");
    }

    // ------------------------------------------------------------------
    // Phase 1: index pass.
    //
    // Build a flat sorted array of (virtual_offset, tag_index) pairs.
    // tag_index is a small integer into a tag-string table, avoiding
    // repeated string storage in the flat array.
    // ------------------------------------------------------------------
    log.step("Phase 1: indexing...");

    samFile*   in  = sam_open(inputBam, "rb");
    if (!in) { std::cerr << "Error: cannot open: " << inputBam << '\n'; return 1; }
    hts_set_threads(in, std::min(nThreads, 4));

    bam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) { std::cerr << "Error: bad header.\n"; sam_close(in); return 1; }

    // Tag string table: string -> uint32_t index.
    std::unordered_map<std::string, uint32_t> tagIndex;
    std::vector<std::string>                  tagStrings; // index -> string

    // Flat index: pairs of (bgzf_virtual_offset, tag_idx).
    // Sorted by offset after the pass so Phase 2 reads sequentially.
    struct OffsetTag { int64_t offset; uint32_t tagIdx; };
    std::vector<OffsetTag> flatIndex;
    flatIndex.reserve(1 << 22); // 4M initial

    uint64_t totalReads    = 0;
    uint64_t assignedReads = 0;
    uint64_t skippedReads  = 0;
    const uint64_t logEvery = 1'000'000;
    uint64_t nextLog = logEvery;

    bam1_t* rec = bam_init1();
    BGZF*   bgzfIn = hts_get_bgzfp(in);

    while (true) {
        int64_t off = bgzf_tell(bgzfIn);
        int r = sam_read1(in, hdr, rec);
        if (r < 0) break;
        ++totalReads;

        if (totalReads >= nextLog) {
            log.step("  indexed " + fmtCount(totalReads)
                   + "  |  tags: " + fmtCount(tagStrings.size()));
            nextLog += logEvery;
        }

        uint8_t* aux = bam_aux_get(rec, splitTag);
        if (!aux) { ++skippedReads; continue; }

        char auxType = bam_aux_type(aux);
        if (auxType != 'Z' && auxType != 'A') { ++skippedReads; continue; }

        const char* tv  = nullptr;
        char singleBuf[2] = {0,0};
        if (auxType == 'Z') {
            tv = bam_aux2Z(aux);
        } else {
            singleBuf[0] = bam_aux2A(aux);
            tv = singleBuf;
        }
        if (!tv || tv[0] == '\0') { ++skippedReads; continue; }

        if (useAllowlist && !allowlist.count(tv)) { ++skippedReads; continue; }

        // Intern the tag string.
        auto [it, inserted] = tagIndex.emplace(tv, static_cast<uint32_t>(tagStrings.size()));
        if (inserted) tagStrings.push_back(tv);

        flatIndex.push_back({off, it->second});
        ++assignedReads;
    }
    bam_destroy1(rec);
    sam_close(in);

    log.step("Phase 1 complete: " + fmtCount(totalReads) + " reads"
           + "  |  assigned: " + fmtCount(assignedReads)
           + "  |  skipped: "  + fmtCount(skippedReads)
           + "  |  unique tags: " + fmtCount(tagStrings.size()));

    // Sort by BGZF virtual offset so Phase 2 reads the input sequentially.
    log.step("Sorting index by file offset...");
    std::sort(flatIndex.begin(), flatIndex.end(),
              [](const OffsetTag& a, const OffsetTag& b){ return a.offset < b.offset; });

    // ------------------------------------------------------------------
    // Phase 2: single sequential scan, write to open output handles.
    //
    // All output files are kept open simultaneously if fd budget allows.
    // If a tag's file handle was evicted due to fd pressure, it is reopened
    // via the raw BGZF-append path (strip EOF, append, write new EOF on close).
    //
    // Because the flat index is sorted by input offset, we walk the input
    // forward-only.  No backward seeks on the input.
    // ------------------------------------------------------------------
    log.step("Phase 2: writing output BAMs (sequential scan)...");

    // Initialise output file state.
    std::vector<OutFile> outFiles(tagStrings.size());
    for (size_t i = 0; i < tagStrings.size(); ++i)
        outFiles[i].path = outPrefix + tagStrings[i] + ".bam";

    // LRU eviction list: indices of currently open output files.
    // We use a simple vector; with thousands of cells, the eviction path
    // is only hit when cell count > fd budget.
    std::vector<uint32_t> openHandles; // indices into outFiles, in open order
    int openCount = 0;

    auto evictOne = [&]() {
        // Close the oldest open handle (front of openHandles).
        if (openHandles.empty()) return;
        uint32_t victim = openHandles.front();
        openHandles.erase(openHandles.begin());
        if (outFiles[victim].fp) {
            sam_close(outFiles[victim].fp);
            outFiles[victim].fp = nullptr;
            --openCount;
        }
    };

    auto ensureOpen = [&](uint32_t idx) -> bool {
        OutFile& of = outFiles[idx];
        if (of.fp) return true; // already open

        // Need to open: evict if at limit.
        while (openCount >= maxOpen) evictOne();

        bool ok = of.headerWritten
                ? reopenForAppend(of)
                : openAndWriteHeader(of, hdr);
        if (!ok) return false;

        openHandles.push_back(idx);
        ++openCount;
        return true;
    };

    // Re-open input for sequential read.
    samFile* in2 = sam_open(inputBam, "rb");
    if (!in2) { std::cerr << "Error: cannot reopen input.\n"; return 1; }
    hts_set_threads(in2, std::min(nThreads, 4));

    bam_hdr_t* hdr2 = sam_hdr_read(in2);
    if (!hdr2) { std::cerr << "Error: bad header on reopen.\n"; sam_close(in2); return 1; }

    bam1_t*   rec2   = bam_init1();
    BGZF*     bgzfIn2 = hts_get_bgzfp(in2);
    int64_t   curOff  = bgzf_tell(bgzfIn2); // position after header
    bool      anyError = false;

    uint64_t written = 0;
    nextLog = logEvery;

    for (const OffsetTag& ot : flatIndex) {
        // Seek only if we are not already at the right position.
        // Because the index is sorted, this seek is almost never needed
        // (consecutive records in same tag may not be adjacent in file,
        //  but we visit them in file order so the next record is always
        //  ahead — the only "seek" is advancing past unneeded reads,
        //  which sam_read1 does naturally).
        if (ot.offset != curOff) {
            if (bgzf_seek(bgzfIn2, ot.offset, SEEK_SET) < 0) {
                std::cerr << "Error: seek failed.\n";
                anyError = true; break;
            }
        }

        if (sam_read1(in2, hdr2, rec2) < 0) {
            std::cerr << "Error: read failed at offset " << ot.offset << '\n';
            anyError = true; break;
        }
        curOff = bgzf_tell(bgzfIn2);

        if (!ensureOpen(ot.tagIdx)) { anyError = true; break; }

        OutFile& of = outFiles[ot.tagIdx];
        if (sam_write1(of.fp, hdr2, rec2) < 0) {
            std::cerr << "Error: write failed to " << of.path << '\n';
            anyError = true; break;
        }
        ++of.nreads;
        ++written;

        if (written >= nextLog) {
            log.step("  written " + fmtCount(written) + " / "
                   + fmtCount(assignedReads) + " reads");
            nextLog += logEvery;
        }
    }

    bam_destroy1(rec2);
    bam_hdr_destroy(hdr2);
    sam_close(in2);

    // Close all open output handles.
    log.step("Closing output files...");
    for (auto& of : outFiles) {
        if (of.fp) { sam_close(of.fp); of.fp = nullptr; }
    }

    // ------------------------------------------------------------------
    // Summary.
    // ------------------------------------------------------------------
    bam_hdr_destroy(hdr);

    uint64_t emptyFiles = 0;
    for (auto& of : outFiles) if (of.nreads == 0) ++emptyFiles;

    log.step("=== Summary ===");
    log.step("  Total reads   : " + fmtCount(totalReads));
    log.step("  Assigned      : " + fmtCount(assignedReads));
    log.step("  Skipped       : " + fmtCount(skippedReads));
    log.step("  Output files  : " + fmtCount(outFiles.size())
           + "  (empty: " + fmtCount(emptyFiles) + ")");
    log.step("  Elapsed       : " + fmtDuration(log.elapsed()));

    if (anyError) { log.step("WARNING: errors occurred."); return 1; }
    return 0;
}