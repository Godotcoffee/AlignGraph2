//
// Created by amadeus on 19-12-30.
//

#include <atomic>
#include <iostream>
#include "kmer/FileKmerIterator.hpp"

int main() {
    FileKmerIterator it("/mnt/hdd1/Project/Genome/counter/yeast.kmer");
    std::atomic_size_t cnt(0);
    std::cout << it.kSize() << std::endl;
    it.iterate([&](std::uint64_t code) {
            ++cnt;
        }, 16);

    std::cout << cnt << std::endl;

    return 0;
}