//
// Created by amadeus on 19-11-12.
//

#include <iostream>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <graph/UnionSet.hpp>
#include <thread>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <mutex>
#include <args.hxx>
#include "align/AlignmentHelper.hpp"
#include "seq/SeqHelper.hpp"


std::unordered_map<std::string, std::set<std::pair<std::string, bool>>>
parseCtgToRef(const std::string &ctgPath, const std::string &ctgToRefPath, std::size_t topK = 2, double ctgToRefTotalRatio = 0.2) {

    std::vector<std::map<std::pair<std::string, bool>, std::vector<bool>>> ctgCover;
    std::unordered_map<std::string, std::size_t> ctgMap;
    std::unordered_map<std::string, std::size_t> refSize;
    std::unordered_map<std::string, std::size_t> ctgSize;

    SeqHelper::autoLoadFromFile(ctgPath, [&](const std::string &name, const std::string &seq) {
        auto shortComment = name.substr(1);
        ctgMap[shortComment] = ctgCover.size();
        ctgSize[shortComment] = seq.size();
        ctgCover.emplace_back();
    });

    AlignmentHelper::simpleLoadFromRefFile(ctgToRefPath, [&](SimpleAlign &align) {
        auto idx = ctgMap[align.queryName];
        auto it = ctgCover[idx].find({align.refName, align.forward});
        if (it == ctgCover[idx].end()) {
            it = ctgCover[idx].emplace(std::make_pair(align.refName, align.forward),
                                       std::vector<bool>(align.querySize, false)).first;
            refSize[align.refName] = align.refSize;
        }

        auto &arr = it->second;
        std::fill(arr.begin() + align.queryBegin, arr.begin() + align.queryEnd, true);

    });

    std::unordered_map<std::string, std::set<std::pair<std::string, bool>>> refWithCtg;

    for (auto &ctg : ctgMap) {
        auto &ctgName = ctg.first;
        auto &cov = ctgCover[ctg.second];

        std::vector<std::pair<std::size_t, std::pair<std::string, bool>>> covToRef;

        for (auto &refCov : cov) {
            auto refIdx = refCov.first;
            auto &record = refCov.second;

            auto total = record.size();
            auto cnt = std::accumulate(record.begin(), record.end(), static_cast<decltype(total)>(0));

            if (cnt * 1.0 / total >= ctgToRefTotalRatio) {
                covToRef.emplace_back(cnt, refIdx);
            }
        }

        std::sort(covToRef.begin(), covToRef.end(), [](const std::pair<std::size_t, std::pair<std::string, bool>> &lhs,
                                                       const std::pair<std::size_t, std::pair<std::string, bool>> &rhs) {
            return lhs.first > rhs.first;
        });

        for (std::size_t i = 0; i < topK && i < covToRef.size(); ++i) {
            refWithCtg[covToRef[i].second.first].insert({ctgName, covToRef[i].second.second});
        }
    }

    decltype(refWithCtg) filterRefWithCtg;

    for (auto &withCtg : refWithCtg) {
        //if ((withCtg.second.size() > 1) || (withCtg.second.size() == 1 && refSize[withCtg.first] >= ctgSize[withCtg.second.begin()->first] * 1.2)) {
        if (withCtg.second.size() > 1) {
            std::unordered_set<std::string> filter;

            for (auto &ctg :withCtg.second) {
                if (filter.count(ctg.first) == 0) {
                    filterRefWithCtg[withCtg.first].insert(ctg);
                    filter.insert(ctg.first);
                }
            }
        }
    }

    return filterRefWithCtg;
}

void filterReadAndRef(const std::string &readPath, const std::string &readOutDir,
                      const std::string &toCtgPath, const std::string &toCtgOutDir,
                      const std::string &toRefPath, const std::string &toRefOutDir,
                      const std::unordered_map<std::string, std::set<std::pair<std::string, bool>>> &filterMap) {
    std::unordered_map<std::string, std::vector<std::size_t>> ctgMap;
    std::unordered_map<std::string, std::vector<std::size_t>> refMap;
    std::map<std::size_t, std::unordered_set<std::size_t>> readMap;
    std::vector<std::set<std::size_t>> readIdx;

    std::vector<std::unique_ptr<std::ofstream>> ctgOfs, refOfs, readOfs;

    for (auto &withCtg : filterMap) {
        refMap[withCtg.first].push_back(refOfs.size());
        for (auto &ctg : withCtg.second) {
            ctgMap[ctg.first].push_back(ctgOfs.size());
        }
        readIdx.emplace_back();
        ctgOfs.emplace_back(new std::ofstream(toCtgOutDir + "/" + std::to_string(ctgOfs.size()) + ".ctg.ref"));
        refOfs.emplace_back(new std::ofstream(toRefOutDir + "/"+ std::to_string(refOfs.size()) + ".ref.ref"));
        readOfs.emplace_back(new std::ofstream(toRefOutDir + "/"+ std::to_string(readOfs.size()) + ".new.fastq"));
    }

    AlignmentHelper::loadFromRefFile(toCtgPath, [&](SimpleAlign &align, const std::string &l1, const std::string &l2, const std::string &l3) {
        auto &arr = ctgMap[align.refName];
        for (auto &idx : arr) {
            std::size_t id = std::stoll(align.queryName);
            readIdx[idx].insert(id);
            readMap[id].insert(idx);

            auto &of = *ctgOfs[idx];
            of << l1 << std::endl;
            of << l2 << std::endl;
            of << l3 << std::endl;
        }
    });

    AlignmentHelper::loadFromRefFile(toRefPath, [&](SimpleAlign &align, const std::string &l1, const std::string &l2, const std::string &l3) {
        auto &arr = refMap[align.refName];
        for (auto &idx : arr) {
            std::size_t id = std::stoll(align.queryName);
            readIdx[idx].insert(id);
            readMap[id].insert(idx);

            auto &of = *refOfs[idx];
            of << l1 << std::endl;
            of << l2 << std::endl;
            of << l3 << std::endl;
        }
    });

    for (auto &of : ctgOfs) {
        of->close();
    }

    for (auto &of : refOfs) {
        of->close();
    }

    for (auto &readSet : readIdx) {
        std::cerr << readSet.size() << std::endl;
    }

    std::size_t count = 0;
    SeqHelper::loadFromFastq(readPath, [&](const std::string &l1, const std::string &l2, const std::string &l3, const std::string &l4) {
        std::size_t id = count + 1;
        auto it = readMap.find(id);
        if (it != readMap.end()) {
            for (auto idx : it->second) {
                auto &of = *readOfs[idx];
                of << "@" << id << std::endl;
                of << l2 << std::endl;
                of << l3 << std::endl;
                of << l4 << std::endl;
            }
        }
        ++count;
        if (count % 50000 == 0) {
            std::cerr << count << std::endl;
        }
    });

    for (auto &of : readOfs) {
        of->close();
    }

    std::cerr << "Total read=" << count << std::endl;
    std::cerr << "Read iterator Done" << std::endl;
}

inline
unsigned acgt(char ch) {
    int n = 0;
    switch (ch) {
        case 'A': case 'a': default:
            n = 0;
            break;
        case 'C': case 'c':
            n = 1;
            break;
        case 'G': case 'g':
            n = 2;
            break;
        case 'T': case 't':
            n = 3;
            break;
    }
    return n;
}

int main(int argc, char *argv[])
{
    args::ArgumentParser parser("pre_process", "");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> readArg(parser, "path", "Read path", {'r', "read"});
    args::ValueFlag<std::string> contigArg(parser, "path", "Contig path", {'c', "contig"});
    args::ValueFlag<std::string> readToCtgArg(parser, "path", "Alignment path of read to contig", {'x', "read_to_ctg"});
    args::ValueFlag<std::string> readToRefArg(parser, "path", "Alignment path of read to reference", {'y', "read_to_ref"});
    args::ValueFlag<std::string> contigToRefArg(parser, "path", "Alignment path of contig to reference", {'z', "ctg_to_ref"});
    args::ValueFlag<std::string> outArg(parser, "path", "Output directory", {'o', "output"});
    args::ValueFlag<std::size_t> topKArg(parser, "unsigned", "Top K of alignment", {'k', "top_k"}, 1);
    args::ValueFlag<double> alignRatioArg(parser, "double", "Alignment min threshold", {'m', "min"}, 0.15);
    args::Flag testArg(parser, "test", "Only do one step", {"test"});

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


    auto ctgTopK                = args::get(topKArg);
    auto ctgToRefTotalRatio     = args::get(alignRatioArg);
    auto readPath               = args::get(readArg);
    auto contigPath             = args::get(contigArg);
    auto readToCtgPath          = args::get(readToCtgArg);
    auto readToRefPath          = args::get(readToRefArg);
    auto ctgToRefPath           = args::get(contigToRefArg);
    auto outDir                 = args::get(outArg);

    auto refWithCtg = parseCtgToRef(contigPath, ctgToRefPath, ctgTopK, ctgToRefTotalRatio);

    for (auto &withCtg : refWithCtg) {
        std::cerr << withCtg.first << std::endl;
        for (auto &ctg : withCtg.second) {
            std::cerr << "\t" << ctg.first << std::endl;
            std::cerr << "\t" << ctg.second << std::endl;
        }
    }

    if (testArg) {
        return 0;
    }

    filterReadAndRef(readPath, outDir,
                     readToCtgPath, outDir,
                     readToRefPath, outDir,
                     refWithCtg);

    std::cerr << "Filter Done" << std::endl;

    std::ofstream of(outDir + "/" + "config.txt");

    std::size_t cnt = 0;
    for (auto &withCtg : refWithCtg) {
        of << withCtg.first << std::endl;
        of << cnt << ".new.fastq" << std::endl;
        of << cnt << ".ctg.ref" << std::endl;
        of << cnt << ".ref.ref" << std::endl;
        for (auto &ctg : withCtg.second) {
            of << ctg.first << std::endl;
            of << ctg.second << std::endl;
        }
        of << std::endl;                    // Empty line
        ++cnt;
    }

    of.close();

    std::cerr << "All Done" << std::endl;

    return 0;
}