#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <htslib/sam.h>
#include <getopt.h>  // For command-line argument parsing
#include <cstring>   // For string manipulation
#include <unordered_map> // For efficient tag-to-reads mapping
#include <fstream>
#include <future>
#include <thread>

std::string addTrailingSlash(const std::string& path) {
    if (path.empty() || path.back() != '/') {
        return path + "/";
    }
    return path;
}

int main(int argc, char* argv[]) {
    const char* usage = "Usage: cell_bam_splitter -i input.bam -o output_dir [-b barcodes.txt] [-t tag]\n";
    const char* input_bam = nullptr;
    const char* output_dir = nullptr;
    const char* barcode_file = nullptr;
    const char* tag = "CB";  // Default tag

    int opt;
    while ((opt = getopt(argc, argv, "i:o:b:t:")) != -1) {
        switch (opt) {
            case 'i':
                input_bam = optarg;
                break;
            case 'o':
                output_dir = optarg;
                break;
            case 'b':
                barcode_file = optarg;
                break;
            case 't':
                tag = optarg;
                break;
            default:
                std::cerr << usage;
                return 1;
        }
    }

    if (!input_bam || !output_dir) {
        std::cerr << usage;
        return 1;
    }

    samFile *in = sam_open(input_bam, "rb");
    if (!in) {
        std::cerr << "Failed to open input BAM file." << std::endl;
        return 1;
    }


    std::set<std::string> tags; // Tag values

    if (barcode_file) {
        std::ifstream barcode_input(barcode_file);
        if (!barcode_input.is_open()) {
            std::cerr << "Failed to open barcode file: " << barcode_file << std::endl;
            return 1;
        }

        std::string barcode;
        while (std::getline(barcode_input, barcode)) {
            tags.insert(barcode);
        }
        barcode_input.close();

    } else { // No barcode_file, gather all unique tags in the file.
        bam_hdr_t *header = sam_hdr_read(in);
        bam1_t *entry = bam_init1();
        size_t read_ct = 0;

        // First pass: Collect unique tag values
        while (sam_read1(in, header, entry) >= 0) {
            ++read_ct;
            uint8_t *tag_data = bam_aux_get(entry, tag);
            if (tag_data) {
                const char* tag_value = bam_aux2Z(tag_data); // Use bam_aux2Z to extract the string value
                tags.insert(tag_value);
            }
        }

        // Close the input BAM file after the first pass
        bam_hdr_destroy(header);
        bam_destroy1(entry);

        // Print the set of unique tags to the screen
        std::cout << "Unique " << tag << " tags found in the input BAM file: " << tags.size() << std::endl;
        std::cout << "Total reads: " << read_ct << std::endl;
    }
    sam_close(in);


    // Map to store the tag-to-reads mapping
    std::unordered_map<std::string, std::vector<bam1_t*>> tagToReadsMap;

    // Collect reads for each tag to be output
    in = sam_open(input_bam, "rb");
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *entry = bam_init1();
    while (sam_read1(in, header, entry) >= 0) {
        uint8_t *tag_data = bam_aux_get(entry, tag);
        if (tag_data) {
            const char* tag_value = bam_aux2Z(tag_data); // Use bam_aux2Z to extract the string value
            if(tags.count(tag_value) > 0) {
                tagToReadsMap[tag_value].push_back(bam_dup1(entry));
            }
        }
    }
    // Close the input BAM file
    sam_close(in);


    size_t totalEntriesInMap = tagToReadsMap.size();
    size_t totalElementsInSet = tags.size();
    size_t elementsIncludedInMap = 0;
    size_t readsIncludedInMap = 0;

    // Iterate over elements in the set and check if they are in the map
    for (const std::string& tag : tags) {
        if (tagToReadsMap.find(tag) != tagToReadsMap.end()) {
            elementsIncludedInMap++;
            readsIncludedInMap += tagToReadsMap[tag].size();
        } else {
           std::cout << "Warning, user-defined " << tag << " was not found in the BAM." << std::endl;
        }
    }

    // Provide a summary to the user
    std::cout << "Summary:" << std::endl;
    std::cout << "Total entries in the map: " << totalEntriesInMap << std::endl;
    std::cout << "Total elements in the set: " << totalElementsInSet << std::endl;
    std::cout << "Elements included in the map: " << elementsIncludedInMap << std::endl;
    std::cout << "Elements not included in the map: " << totalElementsInSet - elementsIncludedInMap << std::endl;
    std::cout << "Total reads in the map: " << readsIncludedInMap << std::endl;

    // Determine the number of CPU threads to use for parallelism
    int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) {
        numThreads = 1; // Fallback to 1 thread if the number of cores can't be determined
    }

    // Create a vector of futures to store the async tasks
    std::vector<std::future<void>> futures;

    // Open the output BAM files and write the reads
    for (const std::string& tag_value : tags) {
        futures.push_back(std::async(std::launch::async, [&]{
            std::string output_bam = addTrailingSlash(output_dir) + tag_value + ".bam";
            std::cout << "Writing " << output_bam << " with " << tagToReadsMap[tag_value].size() << " reads." << std::endl;
            samFile *out = sam_open(output_bam.c_str(), "wb");
            if (!out) {
                std::cerr << "Failed to open output BAM file: " << output_bam << std::endl;
                return;
            }

            if (sam_hdr_write(out, header) < 0) {
                std::cerr << "Failed to write header to " << output_bam << std::endl;
            }

            // Write collected reads for this tag
            const std::vector<bam1_t*>& reads = tagToReadsMap[tag_value];
            for (bam1_t* read : reads) {
                if (sam_write1(out, header, read) < 0) {
                    std::cerr << "Failed to write record to " << tag_value << ".bam" << std::endl;
                }
                bam_destroy1(read); // Clean up the read
            }
            tagToReadsMap[tag_value].clear(); // Free memory
            sam_close(out); // Close the output BAM file
        }));
        if (futures.size() >= numThreads) {
        // Wait for one task to finish before starting a new one to control concurrency
        futures.front().wait();
        futures.erase(futures.begin());
        }
    }

    // Wait for all async tasks to complete
    for (auto& future : futures) {
        future.wait();
    }


    /* for (const std::string& tag_value : tags) {
        std::string output_bam = addTrailingSlash(output_dir) + tag_value + ".bam";
        std::cout << "Writing " << output_bam << " with " << tagToReadsMap[tag_value].size() << " reads." << std::endl;
        samFile *out = sam_open(output_bam.c_str(), "wb");
        if (!out) {
            std::cerr << "Failed to open output BAM file: " << output_bam << std::endl;
            return 1;
        }

        if (sam_hdr_write(out, header) < 0) {
            std::cerr << "Failed to write header to " << output_bam << std::endl;
            return 1;
        }

        // Write collected reads for this tag
        const std::vector<bam1_t*>& reads = tagToReadsMap[tag_value];
        for (bam1_t* read : reads) {
            if (sam_write1(out, header, read) < 0) {
                std::cerr << "Failed to write record to " << tag_value << ".bam" << std::endl;
                return 1;
            }
            bam_destroy1(read); // Clean up the read
        }
        sam_close(out); // Close the output BAM file
    }*/

    // Clean up
    bam_hdr_destroy(header);
    bam_destroy1(entry);

    return 0;
}
