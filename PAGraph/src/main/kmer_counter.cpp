//
// Created by amadeus on 19-12-30.
//

#include <string>
#include <map>
#include <thread>
#include <mutex>
#include <iostream>
#include <fstream>

#include <args.hxx>

#include "seq/SeqHelper.hpp"
#include "kmer/KmerHelper.hpp"
#include "thread/MultiThreadTools.hpp"


void kmerCounter(const std::string &readPath, const std::string &outPath, std::size_t kSize, double threshold,
                 unsigned threadNum = 16, std::size_t hashTablePart = 4, std::size_t partSize = 10240) {
    std::size_t hashTableSize = (1U << (kSize << 1U));

    std::vector<std::vector<std::size_t>> hashTable(hashTablePart, std::vector<std::size_t>(hashTableSize, 0));
    std::vector<std::mutex> mutex(hashTablePart);

    std::vector<std::string> buffer;

    auto action = [&]() {
        std::vector<std::vector<std::uint64_t>> codes(threadNum);

        MultiThreadTools::multiTraversal(buffer, threadNum, [&](unsigned t, std::size_t, std::string &seq) {
            codes[t].clear();
            KmerHelper::kmer2Code(codes[t], seq, kSize);
            auto part = t % hashTablePart;
            {
                std::lock_guard<std::mutex> lockGuard(mutex[part]);
                for (auto &code : codes[t]) {
                    ++hashTable[part][code];
                }
            }
        });

    };

    SeqHelper::autoLoadFromFile(readPath, [&](const std::string &, const std::string &read) {
        buffer.emplace_back(read);
        if (buffer.size() == partSize) {
            action();

            buffer.clear();
        }
    });

    if (!buffer.empty()) {
        action();
    }

    std::map<std::size_t, std::size_t> mergeMap;

    for (std::size_t i = 0; i < hashTableSize; ++i) {
        std::size_t abundance = 0;
        for (std::size_t j = 0; j < hashTablePart; ++j) {
            abundance += hashTable[j][i];
        }
        ++mergeMap[hashTable[0][i] = abundance];
    }

    std::size_t sum = 0;
    std::size_t minAbundance = 0;

    for (auto &kv : mergeMap) {
        sum += kv.second;
        if (1 - sum * 1.0 / hashTableSize <= threshold) {
            minAbundance = kv.first;
            break;
        }
    }

    std::vector<std::vector<std::uint64_t>> availKmer(threadNum);

    MultiThreadTools::multiTraversal(hashTable[0], threadNum, [&](unsigned t, std::size_t i, std::size_t abundance) {
        if (abundance >= minAbundance) {
            availKmer[t].emplace_back(i);
        }
    });

    std::ofstream kmerOf(outPath, std::ios::binary);

    kmerOf.write(reinterpret_cast<const char *>(&kSize), sizeof(std::size_t));

    for (auto &v : availKmer) {
        for (auto &kmer : v) {
            kmerOf.write(reinterpret_cast<const char *>(&kmer), sizeof(std::uint64_t));
        }
    }
}

int main(int argc, char *argv[])
{
    args::ArgumentParser parser("kmer_counter", "");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<unsigned>    threadArg(parser, "thread_num", "Number of thread", {'t', "thread"}, 16);
    args::ValueFlag<std::string> inputArg(parser, "path", "Input of backbone", {'i', "in"});
    args::ValueFlag<std::string> outputArg(parser, "path", "Output of file", {'o', "out"});
    args::ValueFlag<std::size_t> kArg(parser, "k", "kmer size", {'k', "kmer"}, 14);
    args::ValueFlag<double>      thresholdArg(parser, "threshold", "Threshold for min abundance", {'m', "min"}, 0.2);
    args::ValueFlag<std::size_t> partLenArg(parser, "size", "Number of hash table", {'p', "part"}, 4);
    args::ValueFlag<std::size_t> bufferSizeArg(parser, "size", "Size of each batch", {'s', "size"}, 10240);

    if (argc <= 1) {
        std::cerr << parser;
        return 0;
    }

    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cerr << parser;
        return 0;
    } catch (const args::ParseError &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    auto threadNum  = args::get(threadArg);
    auto in         = args::get(inputArg);
    auto out        = args::get(outputArg);
    auto k          = args::get(kArg);
    auto threshold  = args::get(thresholdArg);
    auto partNum    = args::get(partLenArg);
    auto bufferSize = args::get(bufferSizeArg);

    kmerCounter(in, out, k, threshold, threadNum, partNum, bufferSize);

    return 0;
}