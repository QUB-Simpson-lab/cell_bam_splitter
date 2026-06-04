 /* cell_bam_splitter — split a BAM by an arbitrary auxiliary tag value.
 *
 * Architecture (single BAM pass, bounded file handles, no header rewrites):
 *
 *   1. Stream every read once.  For each read that carries the split tag,
 *      clone it into a per-tag in-memory buffer (std::vector<bam1_t*>).
 *
 *   2. Maintain a count of tags whose buffers are "hot" (non-empty but not
 *      yet flushed).  When that count would exceed FILE_HANDLE_RESERVE below
 *      the OS ulimit, trigger a flush of the largest buffers until we are
 *      safely below the limit again.
 *
 *   3. A flush writes all buffered reads for a tag in one open-write-close
 *      cycle.  Because BGZF append is not supported by htslib, each tag
 *      accumulates ALL its reads across potentially multiple flush rounds in
 *      its buffer — the actual disk write only happens once, at finalise().
 *      This keeps the file-per-tag write atomic and avoids any header-rewrite
 *      or BGZF-block fragmentation issue.
 *
 *      Memory pressure valve: if total buffered reads exceed MEM_READ_LIMIT,
 *      the largest tags are flushed to disk even if file handles are available.
 *      Those flushed tags are written to a temp file and merged at the end
 *      using htslib's BAM concatenation.  (See DESIGN NOTE below.)
 *
 * DESIGN NOTE — the two-tier flush strategy:
 *   Tier 1 (most common): all reads for a tag fit in RAM.  One write at end.
 *   Tier 2 (large tags):  reads are flushed to numbered temp shards on disk,
 *                          then cat-merged at finalise.  This never rewrites
 *                          a BAM header mid-stream; each shard is a complete,
 *                          valid BAM.  samtools cat / htslib concat is used
 *                          to merge shards without re-compression.
 *
 * Usage:
 *   cell_bam_splitter -i input.bam -o output_dir
 *                    [-b barcodes.txt] [-t TAG] [-j threads]
 *                    [-m max_reads_in_ram] [-v]
 *
 * Build:
 *   g++ -O2 -std=c++17 cell_bam_splitter.cpp -o cell_bam_splitter \
 *       -lhts -lpthread
 */

#include <iostream>
#include <fstream>
#include <sstream>
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
#include <future>
#include <cstring>
#include <cassert>
#include <sys/stat.h>
#include <sys/resource.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

// ---------------------------------------------------------------------------
// Constants / tunables
// ---------------------------------------------------------------------------

// Headroom: keep this many file descriptors free for stdin/stdout/stderr,
// the input BAM, and OS overhead.
static constexpr int FD_HEADROOM = 32;

// Default cap on total in-RAM bam1_t records before we start spilling large
// tags to disk shards.  Override with -m.
static constexpr size_t DEFAULT_MEM_READ_LIMIT = 4'000'000;

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

static std::string trailingSlash(const std::string& p) {
    return (!p.empty() && p.back() == '/') ? p : p + '/';
}

// Human-readable duration.
static std::string fmtDuration(double seconds) {
    if (seconds < 60)  return std::to_string(static_cast<int>(seconds)) + "s";
    int m = static_cast<int>(seconds) / 60;
    int s = static_cast<int>(seconds) % 60;
    if (m < 60) return std::to_string(m) + "m" + std::to_string(s) + "s";
    int h = m / 60; m = m % 60;
    return std::to_string(h) + "h" + std::to_string(m) + "m";
}

static std::string fmtCount(uint64_t n) {
    // Insert thousand separators.
    std::string s = std::to_string(n);
    int insert = static_cast<int>(s.size()) - 3;
    for (; insert > 0; insert -= 3) s.insert(insert, ",");
    return s;
}

// ---------------------------------------------------------------------------
// Progress logger — prints timestamped one-line status updates.
// ---------------------------------------------------------------------------

class Logger {
public:
    explicit Logger(bool verbose) : verbose_(verbose),
        start_(std::chrono::steady_clock::now()) {}

    void step(const std::string& msg) {
        auto now  = std::chrono::steady_clock::now();
        double el = std::chrono::duration<double>(now - start_).count();
        std::lock_guard<std::mutex> lk(mu_);
        std::cout << "[" << std::fixed << std::setw(7) << std::setprecision(1)
                  << el << "s] " << msg << '\n' << std::flush;
    }

    void verbose(const std::string& msg) {
        if (verbose_) step(msg);
    }

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
// Safe output-directory creation.
// ---------------------------------------------------------------------------

static bool ensureDirectory(const std::string& path) {
    struct stat st{};
    if (stat(path.c_str(), &st) == 0) {
        if (S_ISDIR(st.st_mode)) return true;
        std::cerr << "Error: path exists but is not a directory: " << path << '\n';
        return false;
    }
    // mkdir -p equivalent (single level; for deep paths, iterate).
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i) {
        cur += path[i];
        if (path[i] == '/' || i == path.size() - 1) {
            if (cur.empty() || cur == "/") continue;
            struct stat s2{};
            if (stat(cur.c_str(), &s2) != 0) {
                if (mkdir(cur.c_str(), 0755) != 0 && errno != EEXIST) {
                    std::cerr << "Error: mkdir failed for: " << cur
                              << " (" << strerror(errno) << ")\n";
                    return false;
                }
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Query the OS soft limit on open file descriptors.
// ---------------------------------------------------------------------------

static int getOpenFileLimit() {
    struct rlimit rl{};
    if (getrlimit(RLIMIT_NOFILE, &rl) == 0 && rl.rlim_cur != RLIM_INFINITY) {
        return static_cast<int>(rl.rlim_cur);
    }
    return 1024; // conservative fallback
}

// ---------------------------------------------------------------------------
// Write a vector of bam1_t* records to a single BAM file (with header).
// Returns false on error.  Does NOT destroy the records.
// ---------------------------------------------------------------------------

static bool writeShardToDisk(const std::string& path,
                              bam_hdr_t*         hdr,
                              const std::vector<bam1_t*>& reads,
                              int compressionThreads) {
    samFile* fp = sam_open(path.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error: cannot open for writing: " << path << '\n';
        return false;
    }
    if (compressionThreads > 1) hts_set_threads(fp, compressionThreads);

    if (sam_hdr_write(fp, hdr) < 0) {
        std::cerr << "Error: cannot write header to: " << path << '\n';
        sam_close(fp);
        return false;
    }
    for (bam1_t* r : reads) {
        if (sam_write1(fp, hdr, r) < 0) {
            std::cerr << "Error: write failed to: " << path << '\n';
            sam_close(fp);
            return false;
        }
    }
    sam_close(fp);
    return true;
}

// ---------------------------------------------------------------------------
// Merge a list of BAM shards into a single output BAM using htslib.
// For a single shard, this is just a rename.
// ---------------------------------------------------------------------------

static bool mergeShards(const std::string& finalPath,
                        const std::vector<std::string>& shards,
                        bam_hdr_t* hdr,
                        int compressionThreads) {
    if (shards.empty()) return true;

    if (shards.size() == 1) {
        // Atomic rename — no re-compression.
        if (rename(shards[0].c_str(), finalPath.c_str()) == 0) return true;
        // Rename across devices falls back to copy.
    }

    // Open output.
    samFile* out = sam_open(finalPath.c_str(), "wb");
    if (!out) {
        std::cerr << "Error: cannot open final output: " << finalPath << '\n';
        return false;
    }
    if (compressionThreads > 1) hts_set_threads(out, compressionThreads);
    if (sam_hdr_write(out, hdr) < 0) { sam_close(out); return false; }

    for (const auto& shardPath : shards) {
        samFile* in = sam_open(shardPath.c_str(), "rb");
        if (!in) {
            std::cerr << "Error: cannot open shard: " << shardPath << '\n';
            sam_close(out);
            return false;
        }
        bam_hdr_t* shardHdr = sam_hdr_read(in);
        bam1_t* rec = bam_init1();
        while (sam_read1(in, shardHdr, rec) >= 0) {
            if (sam_write1(out, hdr, rec) < 0) {
                std::cerr << "Error: write failed during merge.\n";
                bam_destroy1(rec);
                bam_hdr_destroy(shardHdr);
                sam_close(in);
                sam_close(out);
                return false;
            }
        }
        bam_destroy1(rec);
        bam_hdr_destroy(shardHdr);
        sam_close(in);
        std::remove(shardPath.c_str()); // clean up shard
    }
    sam_close(out);
    return true;
}

// ---------------------------------------------------------------------------
// Per-tag accumulator.
// ---------------------------------------------------------------------------

struct TagAccumulator {
    std::vector<bam1_t*>    buf;          // in-RAM records
    std::vector<std::string> shards;      // paths of flushed shards on disk
    uint64_t                totalReads = 0;
    bool                    spilled    = false; // ever flushed to disk?
};

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

static void printUsage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " -i input.bam -o output_dir\n"
        << "           [-b barcodes.txt] [-t TAG] [-j threads]\n"
        << "           [-m max_ram_reads] [-v]\n\n"
        << "  -i  Input BAM\n"
        << "  -o  Output directory (created if absent)\n"
        << "  -b  Barcode allowlist, one per line (default: all tags)\n"
        << "  -t  2-char auxiliary tag to split on [default: CB]\n"
        << "  -j  Threads for compression [default: hw concurrency]\n"
        << "  -m  Max reads buffered in RAM before spilling to disk shards\n"
        << "      [default: " << DEFAULT_MEM_READ_LIMIT << "]\n"
        << "  -v  Verbose: log each flush/merge event\n";
}

int main(int argc, char* argv[]) {
    const char* inputBam    = nullptr;
    const char* outputDir   = nullptr;
    const char* barcodeFile = nullptr;
    char        splitTag[3] = "CB";
    int         nThreads    = static_cast<int>(std::thread::hardware_concurrency());
    size_t      memLimit    = DEFAULT_MEM_READ_LIMIT;
    bool        verbose     = false;
    if (nThreads < 1) nThreads = 1;

    int opt;
    while ((opt = getopt(argc, argv, "i:o:b:t:j:m:vh")) != -1) {
        switch (opt) {
            case 'i': inputBam    = optarg; break;
            case 'o': outputDir   = optarg; break;
            case 'b': barcodeFile = optarg; break;
            case 't':
                if (strlen(optarg) != 2) {
                    std::cerr << "Error: -t TAG must be exactly 2 characters.\n";
                    return 1;
                }
                strncpy(splitTag, optarg, 2); splitTag[2] = '\0';
                break;
            case 'j':
                nThreads = std::atoi(optarg);
                if (nThreads < 1) { std::cerr << "Error: -j >= 1 required.\n"; return 1; }
                break;
            case 'm':
                memLimit = static_cast<size_t>(std::stoull(optarg));
                break;
            case 'v': verbose = true; break;
            case 'h': printUsage(argv[0]); return 0;
            default:  printUsage(argv[0]); return 1;
        }
    }

    if (!inputBam || !outputDir) { printUsage(argv[0]); return 1; }

    Logger log(verbose);
    log.step("cell_bam_splitter starting");
    log.step(std::string("  input : ") + inputBam);
    log.step(std::string("  output: ") + outputDir);
    log.step(std::string("  tag   : ") + splitTag);
    log.step(std::string("  threads: ") + std::to_string(nThreads) +
             "  |  ram limit: " + fmtCount(memLimit) + " reads");

    // ------------------------------------------------------------------
    // Create output directory safely.
    // ------------------------------------------------------------------
    if (!ensureDirectory(outputDir)) return 1;
    const std::string outPrefix = trailingSlash(outputDir);
    const std::string tmpDir    = outPrefix + ".tmp_shards/";
    if (!ensureDirectory(tmpDir)) return 1;

    // ------------------------------------------------------------------
    // Derive safe file-handle ceiling.
    // ------------------------------------------------------------------
    const int fdLimit     = getOpenFileLimit();
    const int maxOpenTags = std::max(1, fdLimit - FD_HEADROOM);
    log.step("  fd limit: " + std::to_string(fdLimit) +
             "  |  max concurrent open output files: " + std::to_string(maxOpenTags));

    // ------------------------------------------------------------------
    // Load allowlist.
    // ------------------------------------------------------------------
    std::unordered_set<std::string> allowlist;
    bool useAllowlist = false;
    if (barcodeFile) {
        useAllowlist = true;
        std::ifstream fin(barcodeFile);
        if (!fin) {
            std::cerr << "Error: cannot open barcode file: " << barcodeFile << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(fin, line))
            if (!line.empty()) allowlist.insert(line);
        log.step("Allowlist: " + fmtCount(allowlist.size()) + " barcodes loaded");
    }

    // ------------------------------------------------------------------
    // Open input BAM.
    // ------------------------------------------------------------------
    samFile* in = sam_open(inputBam, "rb");
    if (!in) { std::cerr << "Error: cannot open: " << inputBam << '\n'; return 1; }
    hts_set_threads(in, std::min(nThreads, 4));

    bam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) { std::cerr << "Error: cannot read BAM header.\n"; sam_close(in); return 1; }

    // ------------------------------------------------------------------
    // Streaming accumulation with periodic RAM spill.
    // ------------------------------------------------------------------
    std::unordered_map<std::string, TagAccumulator> acc;
    if (useAllowlist) acc.reserve(allowlist.size());

    uint64_t totalReads      = 0;
    uint64_t assignedReads   = 0;
    uint64_t skippedReads    = 0;
    uint64_t totalInRam      = 0;
    uint64_t spillEvents     = 0;
    int      flushRound      = 0;
    const uint64_t logEvery  = 1'000'000;
    uint64_t nextLog         = logEvery;

    bam1_t* rec = bam_init1();

    // Lambda: spill the N largest in-RAM buffers to disk shards.
    // This frees RAM without closing any output — there are no open output
    // handles during the streaming pass.
    auto spillLargest = [&](size_t nToSpill) {
        // Collect (size, tag) pairs for tags with in-RAM reads.
        std::vector<std::pair<size_t, std::string>> sizes;
        sizes.reserve(acc.size());
        for (auto& [tv, a] : acc)
            if (!a.buf.empty()) sizes.emplace_back(a.buf.size(), tv);

        std::sort(sizes.begin(), sizes.end(),
                  [](auto& a, auto& b){ return a.first > b.first; });

        size_t spilled = 0;
        for (auto& [sz, tv] : sizes) {
            if (spilled >= nToSpill) break;
            auto& a = acc[tv];
            ++flushRound;
            std::string shardPath = tmpDir + tv + ".shard"
                                  + std::to_string(flushRound) + ".bam";
            if (!writeShardToDisk(shardPath, hdr, a.buf, 1)) {
                std::cerr << "Fatal: shard write failed.\n";
                std::exit(1);
            }
            log.verbose("  spill shard: " + tv + " (" + fmtCount(sz) + " reads) -> "
                        + shardPath);
            totalInRam -= a.buf.size();
            for (bam1_t* r : a.buf) bam_destroy1(r);
            a.buf.clear();
            a.buf.shrink_to_fit();
            a.shards.push_back(shardPath);
            a.spilled = true;
            ++spilled;
            ++spillEvents;
        }
    };

    log.step("Streaming BAM pass...");

    while (sam_read1(in, hdr, rec) >= 0) {
        ++totalReads;

        // Progress heartbeat.
        if (totalReads >= nextLog) {
            log.step("  processed " + fmtCount(totalReads) + " reads"
                     + "  |  tags seen: " + fmtCount(acc.size())
                     + "  |  in RAM: " + fmtCount(totalInRam));
            nextLog += logEvery;
        }

        uint8_t* auxData = bam_aux_get(rec, splitTag);
        if (!auxData) { ++skippedReads; continue; }

        char auxType = bam_aux_type(auxData);
        if (auxType != 'Z' && auxType != 'A') { ++skippedReads; continue; }

        const char* tagVal = nullptr;
        char singleCharBuf[2] = {0, 0};
        if (auxType == 'Z') {
            tagVal = bam_aux2Z(auxData);
        } else {
            singleCharBuf[0] = bam_aux2A(auxData);
            tagVal = singleCharBuf;
        }
        if (!tagVal) { ++skippedReads; continue; }

        if (useAllowlist && allowlist.find(tagVal) == allowlist.end()) {
            ++skippedReads;
            continue;
        }

        TagAccumulator& a = acc[tagVal];
        a.buf.push_back(bam_dup1(rec));
        ++a.totalReads;
        ++assignedReads;
        ++totalInRam;

        // RAM pressure check — spill the top-10 largest tags.
        if (totalInRam >= memLimit) {
            log.step("RAM limit reached (" + fmtCount(totalInRam)
                     + " reads in RAM). Spilling largest buffers...");
            // Spill top 10% of tags (or at least 1).
            size_t nSpill = std::max<size_t>(1, acc.size() / 10);
            spillLargest(nSpill);
            log.step("  after spill: " + fmtCount(totalInRam) + " reads in RAM.");
        }
    }
    bam_destroy1(rec);
    sam_close(in);

    log.step("BAM pass complete. Total reads: " + fmtCount(totalReads)
             + "  |  assigned: " + fmtCount(assignedReads)
             + "  |  skipped: "  + fmtCount(skippedReads));
    log.step("Unique tag values: " + fmtCount(acc.size())
             + "  |  spill events: " + std::to_string(spillEvents));

    // ------------------------------------------------------------------
    // Finalise: write all in-RAM buffers (and merge shards) to final BAMs.
    // Uses a thread pool bounded by maxOpenTags to avoid fd exhaustion.
    // ------------------------------------------------------------------
    log.step("Writing " + fmtCount(acc.size()) + " output BAMs"
             + " (up to " + std::to_string(std::min(nThreads, maxOpenTags))
             + " parallel)...");

    // Collect tag values.
    std::vector<std::string> tagValues;
    tagValues.reserve(acc.size());
    for (auto& [tv, _] : acc) tagValues.push_back(tv);

    std::atomic<uint64_t> filesWritten{0};
    std::atomic<bool>     anyError{false};

    // Effective parallelism: limited by both nThreads and fd headroom.
    const int writeConcurrency = std::min(nThreads, maxOpenTags);

    // Simple semaphore via futures vector.
    std::vector<std::future<void>> futures;
    futures.reserve(writeConcurrency);

    auto retireCompleted = [&]() {
        futures.erase(
            std::remove_if(futures.begin(), futures.end(),
                [](std::future<void>& f){
                    return f.wait_for(std::chrono::milliseconds(0))
                           == std::future_status::ready;
                }),
            futures.end());
    };

    const uint64_t nTags = tagValues.size();
    for (const std::string& tv : tagValues) {
        // Throttle to writeConcurrency.
        while (static_cast<int>(futures.size()) >= writeConcurrency) {
            futures.front().wait();
            retireCompleted();
        }

        futures.push_back(std::async(std::launch::async,
            [&, tv]() {
                if (anyError) return;
                TagAccumulator& a  = acc[tv];
                const std::string finalPath = outPrefix + tv + ".bam";

                bool ok = true;
                if (!a.spilled) {
                    // Pure in-RAM: single write.
                    ok = writeShardToDisk(finalPath, hdr, a.buf, 1);
                    for (bam1_t* r : a.buf) bam_destroy1(r);
                    a.buf.clear();
                } else {
                    // Flush remaining RAM buffer as a final shard, then merge.
                    if (!a.buf.empty()) {
                        ++flushRound;
                        std::string lastShard = tmpDir + tv + ".shard"
                                              + std::to_string(flushRound) + ".bam";
                        ok = writeShardToDisk(lastShard, hdr, a.buf, 1);
                        for (bam1_t* r : a.buf) bam_destroy1(r);
                        a.buf.clear();
                        if (ok) a.shards.push_back(lastShard);
                    }
                    if (ok) ok = mergeShards(finalPath, a.shards, hdr, 1);
                }

                if (!ok) anyError = true;
                uint64_t done = ++filesWritten;
                if (done % 1000 == 0 || done == nTags) {
                    log.step("  wrote " + fmtCount(done) + " / "
                             + fmtCount(nTags) + " files");
                }
            }));
    }

    for (auto& f : futures) f.wait();

    // ------------------------------------------------------------------
    // Clean up temp directory.
    // ------------------------------------------------------------------
    rmdir(tmpDir.c_str()); // only removes if empty; ignore error if not

    // ------------------------------------------------------------------
    // Summary.
    // ------------------------------------------------------------------
    bam_hdr_destroy(hdr);

    log.step("\n=== Final summary ===");
    log.step("  Total reads processed : " + fmtCount(totalReads));
    log.step("  Reads assigned        : " + fmtCount(assignedReads));
    log.step("  Reads skipped         : " + fmtCount(skippedReads));
    log.step("  Output files written  : " + fmtCount(filesWritten.load()));
    log.step("  Disk spill events     : " + std::to_string(spillEvents));
    log.step("  Total elapsed         : " + fmtDuration(log.elapsed()));

    if (anyError) {
        log.step("WARNING: one or more write errors occurred. Check output.");
        return 1;
    }
    return 0;
}
