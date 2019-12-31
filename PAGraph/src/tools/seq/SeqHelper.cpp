//
// Created by amadeus on 19-11-13.
//

#include "SeqHelper.hpp"
#include <fstream>

void SeqHelper::loadFromFastq(const std::string &fastqPath,
                              std::function<void(const std::string &, const std::string &, const std::string &,
                                                 const std::string &)> f) {
    std::ifstream inFastq(fastqPath);
    if (inFastq.is_open()) {
        std::string line1, line2, line3, line4;

        for (std::size_t lineCount = 0; ; ++lineCount) {
            if (lineCount % 4 == 0) {
                if (!std::getline(inFastq, line1)) { break; }
            } else if (lineCount % 4 == 1) {
                if (!std::getline(inFastq, line2)) { break; }
            } else if (lineCount % 4 == 2) {
                if (!std::getline(inFastq, line3)) { break; }
            } else {
                if (!std::getline(inFastq, line4)) { break; }
                f(line1, line2, line3, line4);
            }
        }
        inFastq.close();
    }
}

void SeqHelper::loadFromFasta(const std::string &fastaPath,
        std::function<void(const std::string &, const std::string &)> f) {
    std::ifstream inFasta(fastaPath);
    if (inFasta.is_open()) {
        std::string line;
        std::string name;
        std::string buffer;
        while (std::getline(inFasta, line)) {
            if (!line.empty() && line[0] == '>') {
                if (!name.empty()) {
                    f(name, buffer);
                    buffer.clear();
                }
                name = line;
            } else {
                buffer += line;
            }
        }

        if (!name.empty()) {
            f(name, buffer);
            buffer.clear();
        }

        inFasta.close();
    }
}

void SeqHelper::autoLoadFromFile(
        const std::string &path,
        std::function<void(const std::string &,const std::string &)> f) {
    auto fileType = testFileType(path);
    bool useFasta = fileType == "fasta";

    if (useFasta) {
        loadFromFasta(path, [&](const std::string &name, const std::string &seq) {
            f(name, seq);
        });
    } else {
        loadFromFastq(path, [&](const std::string &l1, const std::string &l2, const std::string &l3, const std::string &l4) {
            f(l1, l2);
        });
    }
}

std::string SeqHelper::testFileType(const std::string &path) {
    std::ifstream in(path);

    if (!in) {
        return "unknown";
    }

    std::string firstLine;

    if (!std::getline(in, firstLine) || firstLine.empty()) {
        return "unknown";
    }

    auto first = firstLine.front();

    switch (first) {
        case '@':
            return "fastq";
        case '>': case ';':
            return "fasta";
        default:
            return "unknown";
    }
}
