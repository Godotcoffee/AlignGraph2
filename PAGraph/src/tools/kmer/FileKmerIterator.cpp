#include <utility>
#include <fstream>
#include <thread/MultiThreadTools.hpp>

//
// Created by amadeus on 19-12-26.
//

#include "FileKmerIterator.hpp"

FileKmerIterator::FileKmerIterator(std::string path) : _k(0),  _filePath(std::move(path)) {
    std::ifstream in(_filePath, std::ios::binary);
    in.read(reinterpret_cast<char *>(&_k), sizeof(std::size_t));
}

void FileKmerIterator::iterate(std::function<void(std::uint64_t)> f, unsigned threadNum) const {
    std::ifstream in(_filePath, std::ios::binary);

    std::string line;
    std::vector<std::uint64_t > lines;
    std::uint64_t buffer = 0;

    std::size_t maxBufferSize = 10240;

    auto action = [&]() {
        MultiThreadTools::multiTraversal(lines, threadNum, [&](unsigned, std::size_t, std::uint64_t kmer) {
            f(kmer);
        });
    };

    while (in.read(reinterpret_cast<char *>(&buffer), sizeof(std::uint64_t))) {
        lines.emplace_back(buffer);

        if (lines.size() >= maxBufferSize) {
            action();

            lines.clear();
        }
    }

    if (!lines.empty()) {
        action();
    }
}
