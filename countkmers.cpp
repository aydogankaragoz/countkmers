#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <algorithm>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <string>

#include "3rdParty/bloom_filter.hpp"
#include "3rdParty/google/sparse_hash_map"

void handle_error(const char* msg) {
    perror(msg);
    exit(255);
}

const char* map_file(const char* fname, size_t& length) {
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    return addr;
}

struct eqstr {
    bool operator()(std::string s1, std::string s2) const {
        return (s1 == s2) || ((s1.compare(s2)) == 0);
    }
};

bloom_filter build_bloom_filter(unsigned long long int projected_element_count,
    double false_positive_probability) {

        bloom_parameters parameters;
        parameters.projected_element_count = projected_element_count;
        parameters.false_positive_probability = false_positive_probability;
        parameters.random_seed = 0xA5A5A5A5;
        parameters.compute_optimal_parameters();

        return bloom_filter(parameters);
}

template<typename T>
void printHead(T s, int c) {
    std::cout << "inserting into vector" << "\n";
    std::vector<std::pair<std::string, int> > items;
    for (auto const& kmer : s) {
        items.push_back(std::make_pair(kmer.first, kmer.second));
    }

    // sorting the vector
    std::cout << "sorting the vector" << "\n";
    std::sort(items.begin(), items.end(),
        [](const std::pair<std::string, int>& firstElem, const std::pair<std::string, int>& secondElem) -> bool
        {
            return firstElem.second > secondElem.second;
        });

    for (int i=0; i < c; i++)
        std::cout << items[i].second << " -> " << items[i].first << std::endl;
}
int main(int argc, char *argv[] ) {
    if (argc < 4) {
       std::cerr << "USAGE: ./countkmers <filename> <kmersize> <topcount>\n";
       return -1;
    }

    char* file_name = argv[1];
    int kmersize = std::atoi(argv[2]);
    int topcount = std::atoi(argv[3]);

    char *kmerArray = new char[kmersize];
    int kacKmer = 0;

    bloom_filter filter = build_bloom_filter(100000000, 0.0001);

    google::sparse_hash_map<std::string, int, std::hash<std::string>, eqstr> kmers;

    size_t length;
    auto f = map_file(file_name, length);
    auto l = f + length;
    auto previous_f = f;

    uintmax_t m_numLines = 0;
    while (f && f != l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))) {
        m_numLines++;

        //We are interested in every 2nd line of the 4 line set
        if (m_numLines%4 == 2) {
            auto sequence_len = f - previous_f;

            for (int i = 0; i <= sequence_len - kmersize; i++) {
                //Directly copy k-mer from file
                memcpy(kmerArray, previous_f + i, kmersize);
                kacKmer++;

                std::string kmer_str(kmerArray);
                if (!filter.contains(kmer_str)) {
                    filter.insert(kmer_str);
                } else {
                    kmers[kmer_str]++;
                }
            }
        }
        f++;
        previous_f = f;
    }

    delete [] kmerArray;

    printHead(kmers, topcount);
}
