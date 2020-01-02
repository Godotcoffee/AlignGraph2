#include <iostream>
#include <thread>
#include <atomic>
#include <fstream>
#include <args.hxx>
#include "seq/AutoSeqDatabase.hpp"
#include "thread/MultiThreadTools.hpp"
#include "cns/Alignment.hpp"
#include "cns/AlignData.hpp"
#include "cns/AlnGraphBoost.hpp"

int main(int argc, char *argv[]) {
    /*cxxopts::Options options("pbdgatest", "");
    options.add_options()
            ("t,thread", "Number of thread", cxxopts::value<unsigned>()->default_value("8"))
            ("l,len", "Length of part", cxxopts::value<std::size_t>()->default_value("500"))
            ("i,in", "Input of backbone", cxxopts::value<std::string>())
            ("o,out", "Output of file", cxxopts::value<std::string>())
            ("a,align", "Alignments file", cxxopts::value<std::string>());
    auto optResult = options.parse(argc, argv);*/


    args::ArgumentParser parser("pa_cns", "");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<unsigned> threadArg(parser, "thread_num", "number of thread", {'t', "thread"}, 16);
    args::ValueFlag<std::size_t> partLenArg(parser, "size", "length of part", {'l', "len"}, 10000);
    args::ValueFlag<std::size_t> topKArg(parser, "k", "length of part", {'k', "top_k"}, 3000);
    args::ValueFlag<std::string> inputArg(parser, "path", "input of backbone", {'i', "in"});
    args::ValueFlag<std::string> outputArg(parser, "path", "output of file", {'o', "out"});
    args::ValueFlag<std::string> alignArg(parser, "path", "alignments file", {'a', "align"});

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
    auto partLen    = args::get(partLenArg);
    auto align      = args::get(alignArg);
    auto topK       = args::get(topKArg);

    AutoSeqDatabase seqDatabase(in);

    if (seqDatabase.size() == 0) {
        std::cerr << "[Error] No backbone" << std::endl;
        return 2;
    }

    auto name = seqDatabase[0].name();
    auto backbone = seqDatabase[0].toString();

    std::size_t partNum = (backbone.size() + partLen - 1) / partLen;

    std::cout << "PartNum=" << partNum << std::endl;

    std::vector<std::string> consensusResults(partNum);

    auto alignments = AlignData::readFromRefFile(align, backbone, partLen);

    auto sortAndTopFun = [&](std::size_t i) {
        std::sort(alignments[i].begin(), alignments[i].end(),
                  [](const std::pair<dagcon::Alignment, unsigned> &lhs,
                     const std::pair<dagcon::Alignment, unsigned> &rhs) {
                      return lhs.second > rhs.second;
                  });
        if (alignments[i].size() > topK) {
            alignments[i].resize(topK);
        }
    };

    auto preThreadFun = [&](unsigned startThread) {
        for (std::size_t i = startThread; i < partNum; i += threadNum) {
            sortAndTopFun(i);
        }
    };

    MultiThreadTools::mulitFunctor(threadNum, preThreadFun);

    std::atomic_size_t consensusLen(0);

    auto consensusFun = [&](std::size_t i) {
        //std::cout << i << std::endl;
        //std::cout << "part " << i << "/" << partNum << std::endl;
        auto left = i * partLen;
        auto right = std::min((i + 1) * partLen, backbone.size());

        auto partSkeleton = backbone.substr(left, right - left);

        std::size_t maxCov = partLen;

        //std::cout << "\t" << alignments[i].size() << std::endl;

        auto weights = AlignData::weightAln(alignments[i]);

        auto g = AlnGraphBoost(partSkeleton);
        for (std::size_t ii = 0; ii < weights.size(); ++ii) {
            auto countDown = weights[ii];
            /*while (countDown--) {
                g.addAln(alignments[i][ii].first);
            }*/
            g.addAln(alignments[i][ii].first, countDown);
        }

        g.mergeNodes();

        consensusResults[i] = g.consensus();

        consensusLen += consensusResults[i].size();
    };

    auto threadFun = [&](unsigned startThread) {
        for (std::size_t i = startThread; i < partNum; i += threadNum) {
            consensusFun(i);
        }
    };

    std::vector<std::thread> threads;

    for (unsigned t = 0; t < threadNum; ++t) {
        threads.emplace_back(threadFun, t);
    }

    for (auto &thread : threads) {
        thread.join();
    }

    std::cout << consensusLen << std::endl;
    std::cout << backbone.size() << std::endl;

    const std::size_t lineSize = 70;

    std::ofstream of(out);

    of << ">" << name << std::endl;

    std::size_t cnt = 0;

    for (auto &seq : consensusResults) {
        for (auto ch : seq) {
            of << ch;
            if (++cnt % lineSize == 0) {
                of << std::endl;
                cnt = 0;
            }
        }
    }
    if (cnt > 0) {
        of << std::endl;
    }

    return 0;
}