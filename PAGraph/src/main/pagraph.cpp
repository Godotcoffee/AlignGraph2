#include <unordered_map>
#include <map>
#include <fstream>
#include <regex>
#include <chrono>
#include <utility>
#include <args.hxx>
#include <align/MummerAlignDatabaseV2.hpp>
#include <graph/UnionSet.hpp>
#include "seq/AutoSeqDatabase.hpp"
#include "graph/PABruijnGraph.hpp"
#include "graph/PAlgorithm.hpp"
#include "align/Aligner.hpp"
#include "position/PositionProcessor.hpp"
#include "tools/MyTools.hpp"
#include "align/AlignReference.hpp"
#include "node/KMerAdjNode.hpp"
#include "graph/PAssembly.hpp"
#include "kmer/FileKmerIterator.hpp"

struct Config {
    std::string ref;
    std::vector<std::pair<std::string, bool>> contigs;
    std::string readPath;
    std::string ctgAlnPath;
    std::string refAlnPath;
};

std::vector<Config> loadFromConfig(const std::string &path) {
    std::vector<Config> configs;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line)) {
        Config cfg;
        cfg.ref = line;
        std::getline(in, cfg.readPath);
        std::getline(in, cfg.ctgAlnPath);
        std::getline(in, cfg.refAlnPath);

        while (std::getline(in, line) && !line.empty()) {
            cfg.contigs.emplace_back(line, false);
            std::getline(in, line);
            std::stringstream(line) >> cfg.contigs.back().second;
        }

        configs.push_back(cfg);
    }
    return configs;
}

std::ostream &operator<<(std::ostream &os, const Config &config) {
    os << config.ref << " " << std::endl;
    for (auto &ctg : config.contigs) {
        os << ctg.first << "(" << (ctg.second ? "forward" : "reverse") << ")" << " | ";
    }
    os << std::endl;
    os << config.readPath << std::endl;
    os << config.ctgAlnPath << std::endl;
    os << config.refAlnPath << std::endl;
    return os;
}

void printConfigs(const std::vector<Config> &configs) {
    for (auto &config : configs) {
        std::cout << config << std::endl;
    }
}

int run2(int argc, char* argv[])
{
    args::ArgumentParser parser("pagraph", "");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<unsigned> threadArg(parser, "thread_num", "number of thread", {'t', "thread"}, 16);
    args::ValueFlag<std::string> kmerArg(parser, "path", "solid kmer set path", {'k', "kmer"});
    args::ValueFlag<std::string> readArg(parser, "path", "read path", {'r', "read"});
    args::ValueFlag<std::string> ctgArg(parser, "path", "contig path", {'c', "contig"});
    args::ValueFlag<std::string> refArg(parser, "path", "reference path", {'R', "ref"});
    args::ValueFlag<std::string> preArg(parser, "path", "pre process directory", {'p', "pre_process"});
    args::ValueFlag<std::string> ctgToRefArg(parser, "path", "alignment path of contig to reference", {'a', "aln"});
    args::ValueFlag<std::string> outArg(parser, "path", "output directory", {'o', "output"});
    args::ValueFlag<std::size_t> minLenArg(parser, "len", "minimum path length", {'l', "length"}, 50);
    args::ValueFlag<std::size_t> devArg(parser, "dist", "distance to join vertices", {"epsilon"}, 10);
    args::ValueFlag<std::size_t> covArg(parser, "cov", "coverage to filter", {'v'}, 1);

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

    auto threadNum = args::get(threadArg);
    auto kmerPath = args::get(kmerArg);
    auto readsPath = args::get(readArg);
    auto contigsPath = args::get(ctgArg);
    auto referencesPath = args::get(refArg);
    auto inputDir = args::get(preArg);
    auto contigToRefPath = args::get(ctgToRefArg);
    auto outDir = args::get(outArg);

    int readToCtgTopK = -1;
    int ctgToRefTopK = -1;
    int readToRefTopK = -1;
    int outerSample = 3;
    int innerSample = 1;
    double readToCtgRatio = 0.35;
    double readToRefRatio = 0.10;
    double ctgToRefRatio = 0.00;
    double ctgToRefTotalRatio = 0.1;
    std::size_t ctgToRefMinLen = 50;
    //std::size_t posError = 10;
    std::size_t posError = args::get(devArg);
    std::size_t covFilter = args::get(covArg);
    double errorRate = 0.15;
    double startSplit = 0.90;
    auto minLen = args::get(minLenArg);

    std::cout << "Loading config.txt" << std::endl;

    auto configs = loadFromConfig(inputDir + "/config.txt");
    printConfigs(configs);

    std::cout << "Loading KMer" << std::endl;

    auto pKmerIt = std::make_shared<FileKmerIterator>(kmerPath);

    std::cout << "Done! k=" << pKmerIt->kSize() << std::endl;

    std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

    std::cout << "Loading Contigs" << std::endl;

    auto pContigDB = std::make_shared<AutoSeqDatabase>(contigsPath);

    std::cout << "Done! contigs number=" << pContigDB->size() << std::endl;

    std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

    std::cout << "Loading References" << std::endl;

    auto pRefDB = std::make_shared<AutoSeqDatabase>(referencesPath);

    AutoSeqDatabase as = *pRefDB;

    std::cout << "Done! reference number=" << pRefDB->size() << std::endl;

    std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

    std::cout << "Loading ContigToRef" << std::endl;

    auto pContigToRef = std::make_shared<MummerAlignDatabaseV2>(contigToRefPath);

    std::cout << "Done! number=" << pContigToRef->size() << std::endl;

    std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

    std::cout << "Building original pa Graph" << std::endl;

    auto pPaGraph = std::make_shared<PABruijnGraph>(*pKmerIt, threadNum);

    std::cout << "Done! kmer number=" << pPaGraph->availableKmerNumber() << std::endl;

    std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

    std::cout << "Done!" << std::endl;

    std::unordered_set<std::string> doneCtg;
    std::unordered_set<std::string> okCtg;

    std::size_t helpCnt = 0;

    for (auto &config : configs) {
        pPaGraph->resetAllNodes(threadNum);
        std::cout << "Use Ref: " << config.ref << std::endl;
        std::cout << "Loading read from " << config.readPath << std::endl;

        auto pReadDB = std::make_shared<AutoSeqDatabase>(inputDir + "/" + config.readPath);

        std::cout << "Done! reads number=" << pReadDB->size() << std::endl;
        std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

        std::cout << "Loading read to contig from " << config.ctgAlnPath << std::endl;

        auto pReadToContig = std::make_shared<MecatAlignDatabase>(inputDir + "/" + config.ctgAlnPath);

        std::cout << "Done! aln number=" << pReadToContig->size() << std::endl;
        std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

        std::cout << "Loading read to ref from " << config.refAlnPath << std::endl;

        auto pReadToRef = std::make_shared<MecatAlignDatabase>(inputDir + "/" + config.refAlnPath);

        std::cout << "Done! aln number=" << pReadToRef->size() << std::endl;
        std::cout << "Mem=" << MyTools::totalVirtualMem() << std::endl;

        PositionProcessor pp(pPaGraph, pReadDB, pContigDB, pRefDB, pReadToContig, pReadToRef, pContigToRef);

        pp.setReadToCtgTopK(readToCtgTopK);
        pp.setReadToRefTopK(readToRefTopK);
        pp.setCtgToRefTopK(ctgToRefTopK);
        pp.setOuterSample(outerSample);
        pp.setInnerSample(innerSample);
        pp.setPositionError(posError);
        pp.setReadToCtgRatio(readToCtgRatio);
        pp.setReadToRefRatio(readToRefRatio);
        pp.setCtgToRefRatio(ctgToRefRatio);
        pp.setCtgToRefTotalRatio(ctgToRefTotalRatio);
        pp.setCtgToRefMinLen(ctgToRefMinLen);
        pp.setCovFilter(covFilter);
        pp.setThreadNum(threadNum);

        pp.clearRefFilter(false);
        pp.clearCtgFilter(false);

        pp.setRefFilter(config.ref, true);

        std::set<std::pair<std::string, bool>> usedCtg;

        for (auto &ctg : config.contigs) {
            pp.setCtgFilter(ctg.first, ctg.second, true);
            usedCtg.emplace(ctg);
        }

        std::cout << "Pre Process" << std::endl;

        pp.preProcess();

        std::cout << "Process Mem=" << MyTools::totalVirtualMem() << std::endl;

        pp.process();

        auto successCtg = PAssembly::testTravel5(
                outDir,
                std::to_string(helpCnt) + "_",
                pPaGraph,
                pReadDB,
                pContigDB,
                pRefDB,
                std::make_shared<PositionMapper>(*pContigDB),
                std::make_shared<PositionMapper>(*pRefDB),
                usedCtg,
                posError * 2,
                errorRate,
                startSplit,
                minLen,
                threadNum
        );

        ++helpCnt;

        for (auto &success : successCtg) {
            okCtg.emplace(success.first);
        }
    }

    std::ofstream ctgList(outDir + "/contig.txt");

    for (auto &processCtg : okCtg) {
        ctgList << processCtg << std::endl;
    }

    return EXIT_SUCCESS;
}

std::ostream &operator<<(std::ostream &os, const AlignInf &align) {
    os << align.getQueryName() << "\t" << align.getRefName() << "\t"
    << align.getQueryBegin() << "\t" << align.getQueryEnd() << "\t"
    << align.getRefBegin() << "\t" << align.getRefEnd() << "\t" << (align.isForward() ? "F" : "R");
    return os;
}

template<typename T>
void combatSeq(const std::vector<std::vector<std::size_t>> &seqs,
               std::size_t start, bool forward, T functor) {
    std::size_t nextIdx = start;
    std::size_t nextPos = 0;

    std::set<std::size_t> ctgSet;
    ctgSet.emplace(nextIdx);

    while (true) {
        bool go = functor(nextIdx);
        if (!go) {
            break;
        }

        if (seqs[nextIdx].empty()) {
            break;
        }

        nextIdx = seqs[nextIdx].back();

        if (ctgSet.count(nextIdx) > 0) {
            break;
        }

        ctgSet.emplace(nextIdx);
    }
}

void test1() {
    std::set<std::pair<std::string, bool>> ctgSet{
            {"A", true},
            {"B", true},
            {"C", true},
            {"D", true},
            {"E", true},
            {"F", true},
            {"G", true}
    };

    std::map<std::string, std::size_t> ctgMap {
            {"A", 0},
            {"B", 1},
            {"C", 2},
            {"D", 3},
            {"E", 4},
            {"F", 5},
            {"G", 6}
    };

    std::vector<std::pair<std::string, std::size_t>> t {
            {"A", 130},
            {"B", 50},
            {"C", 150},
            {"D", 10},
            {"E", 10},
            {"F", 10},
            {"G", 100}
    };

    std::vector<std::size_t> inDegrees {
            1, 2, 0, 1, 0, 1, 0
    };

    std::vector<std::vector<std::size_t>> cd {
            {1}, {3}, {1}, {}, {5}, {}, {0}
    };

    std::cout << "[Union] Start" << std::endl;

    std::map<std::pair<std::string, bool>, std::size_t> ctgIdxHelper;
    std::vector<std::pair<std::string, bool>> ctgTable;

    for (auto &ctgName : ctgSet) {
        ctgIdxHelper[ctgName] = ctgTable.size();
        ctgTable.emplace_back(ctgName);
    }

    std::vector<bool> touched(ctgIdxHelper.size(), false);
    UnionSet unionSet;
    unionSet.init(ctgIdxHelper.size());

    for (auto &ctgName : ctgSet) {
        auto ctgIdx = ctgMap[ctgName.first];
        auto i = ctgIdx;
        if (inDegrees[i] > 0) { continue; }

        auto mainIdx = ctgIdxHelper[ctgName];
        touched[mainIdx] = true;

        combatSeq(cd, i, true,
                  [&](std::size_t ctgId) -> bool {
                      if (ctgId != ctgIdx) {
                          auto helperIdx = ctgIdxHelper[{t[ctgId].first, true}];
                          unionSet.unionTwo(helperIdx, mainIdx);
                          if (touched[helperIdx]) {
                              return false;
                          } else {
                              touched[helperIdx] = true;
                              return true;
                          }
                      } else {
                          return true;
                      }
                  });
    }

    std::vector<std::vector<std::size_t>> unionMerge(ctgIdxHelper.size());
    for (std::size_t i = 0; i < ctgTable.size(); ++i) {
        unionMerge[unionSet.find(i)].emplace_back(i);
    }

    std::set<std::pair<std::string, bool>> filterCtgSet;

    for (auto &eachSet : unionMerge) {
        if (eachSet.empty()) { continue; }
        std::size_t maxSize = 0;
        std::size_t chosenOne = eachSet.front();
        for (auto &idx : eachSet) {
            auto &ctgName = ctgTable[idx];
            auto ctgIdx = ctgMap[ctgName.first];
            auto i = ctgIdx;
            if (inDegrees[i] > 0) { continue; }

            std::size_t len = 0;

            combatSeq(cd, i, true,
                      [&](std::size_t ctgId) -> bool {
                          len += t[ctgId].second;
                          return true;
                      });

            if (len > maxSize) {
                maxSize = len;
                chosenOne = idx;
            }
        }
        filterCtgSet.emplace(ctgTable[chosenOne]);
    }


    std::cout << "[Union] End" << std::endl;

    for (auto &ctgSet : filterCtgSet) {
        std::cout << ctgSet.first << std::endl;
    }
}


int main (int argc, char* argv[])
{
    return run2(argc, argv);
}
