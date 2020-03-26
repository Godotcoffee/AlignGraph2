//
// Created by amadeus on 19-5-13.
//

#include "PAssembly.hpp"
#include "UnionSet.hpp"

#include <map>

std::set<std::pair<std::string, bool>>
PAssembly::testTravel5(const std::string &outDir, const std::string &prefix, std::shared_ptr<PABruijnGraph> pGraph,
                       std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
                       std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
                       std::shared_ptr<const ISeqDatabase<SeqInf>> pRefDB,
                       std::shared_ptr<const PositionMapper> pCtgMapper,
                       std::shared_ptr<const PositionMapper> pRefMapper,
                       const std::set<std::pair<std::string, bool>> &ctgSet, std::size_t deviation, double errorRate,
                       double startSplit, unsigned threadNum) {

    std::set<std::pair<std::string, bool>> successCtgSet;

    PAlgorithm algo(pGraph, pReadDB, pContigDB, pRefDB, pCtgMapper, pRefMapper, threadNum);

    std::vector<PAlgorithm::TravelSequence> results(pContigDB->size() * 2);
    std::vector<std::size_t> inDegrees(pContigDB->size() * 2);

    for (auto &ctgName : ctgSet) {
        //for (std::size_t i = 0; i < pContigDB->size(); ++i) {
        auto ctgIdx = pContigDB->seqId(ctgName.first);
        auto ctgOffset = ctgName.second ? 0 : 1;

        std::cout << "[Travel] " << ctgIdx << " - " << (*pContigDB)[ctgIdx].name() << " - " << (*pContigDB)[ctgIdx].size() << std::endl;
        //for (int j = 0; j < 2; ++j) {
        std::cout << "[Travel] " << (ctgOffset == 0 ? "forward" : "reverse") << std::endl;
        results[2 * ctgIdx + ctgOffset] = algo.travelSequence(ctgIdx, ctgOffset == 0, deviation, errorRate, startSplit);

        std::ofstream of(outDir + "/" + prefix
                         + std::to_string(ctgIdx) + "_" + std::to_string(ctgOffset) + ".txt");

        of << ctgName.first << "\t" << (*pContigDB)[ctgIdx].size() << std::endl;

        for (auto &s : results[2 * ctgIdx + ctgOffset]) {
            auto firstPos = s.first.getPosition().first;
            auto secondPos = s.first.getPosition().second;
            auto dualPos1 = pCtgMapper->singleToDual(firstPos);
            auto dualPos2 = pRefMapper->singleToDual(secondPos);

            of << pGraph->toString(s.first) << "\t" << s.second
               << "\t" << dualPos1.first << "," << dualPos1.second
               << "\t" << dualPos2.first << "," << dualPos2.second << std::endl;
        }

        if (PAlgorithm::seqSize(results[2 * ctgIdx + ctgOffset]) < (*pContigDB)[ctgIdx].size() * startSplit * 0.9) {
            results[2 * ctgIdx + ctgOffset].clear();
        }

        if (!results[2 * ctgIdx + ctgOffset].empty()) {
            auto lastCtgPos = results[2 * ctgIdx + ctgOffset].back().first.getPosition().first;

            if (lastCtgPos != 0) {
                auto dualPos = pCtgMapper->singleToDual(lastCtgPos);
                auto idx = static_cast<decltype(ctgIdx)>(std::abs(dualPos.first) - 1);
                auto forward = dualPos.first > 0 ? 0 : 1;

                if ((idx != ctgIdx || forward != ctgOffset)) {
                    ++inDegrees[2 * idx + forward];
                }
            }
        }

        //}
        std::cout << "[Travel] End" << std::endl;
    }

    for (auto &ctgName : ctgSet) {
        auto ctgIdx = pContigDB->seqId(ctgName.first);
        auto ctgOffset = ctgName.second ? 0 : 1;

        if (!results[2 * ctgIdx + ctgOffset].empty()) {
            auto lastCtgPos = results[2 * ctgIdx + ctgOffset].back().first.getPosition().first;

            if (lastCtgPos != 0) {
                auto dualPos = pCtgMapper->singleToDual(lastCtgPos);
                auto idx = static_cast<decltype(ctgIdx)>(std::abs(dualPos.first) - 1);
                auto forward = dualPos.first > 0 ? 0 : 1;

                if ((idx != ctgIdx || forward != ctgOffset)) {
                    if (results[2 * idx + forward].empty()) {
                        results[2 * ctgIdx + ctgOffset].pop_back();
                        --inDegrees[2 * idx + forward];
                    }
                }
            }
        }
    }

    for (std::size_t i = 0; i < inDegrees.size(); ++i) {
        if (inDegrees[i] > 0) {
            std::cout << "\t" << (i / 2) << " " << (i % 2 == 0 ? "forward" : "reverse") << " " << inDegrees[i]
                      << std::endl;
        }
    }

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
        auto ctgIdx = pContigDB->seqId(ctgName.first);
        auto i = ctgIdx * 2 + (ctgName.second ? 0 : 1);
        if (inDegrees[i] > 0 || results[i].empty()) { continue; }

        auto mainIdx = ctgIdxHelper[ctgName];
        touched[mainIdx] = true;

        combatSeq(results, *pCtgMapper, i / 2, i % 2 == 0,
                  [&](std::size_t ctgId, bool forward, std::size_t startPos) -> bool {
                      if (ctgId != ctgIdx || forward != ctgName.second) {
                          auto helperIdx = ctgIdxHelper[{(*pContigDB)[ctgId].name(), forward}];
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
            auto ctgIdx = pContigDB->seqId(ctgName.first);
            auto i = ctgIdx * 2 + (ctgName.second ? 0 : 1);
            if (inDegrees[i] > 0 || results[i].empty()) { continue; }

            std::size_t len = 0;

            combatSeq(results, *pCtgMapper, i / 2, i % 2 == 0,
                      [&](std::size_t ctgId, bool forward, std::size_t startPos) -> bool {
                          len += (*pContigDB)[ctgId].size();
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
    
    std::cout << "Start From:" << std::endl;
    for (auto &ctgName : filterCtgSet) {
        std::cout << "\t" << ctgName.first << " " << ctgName.second << std::endl;
    }

    std::cout << "[Assembly] Start" << std::endl;

    std::size_t nameCnt = 0;

    for (auto &ctgName : filterCtgSet) {
        //for (std::size_t i = 0; i < results.size(); ++i) {
        auto ctgIdx = pContigDB->seqId(ctgName.first);
        auto i = ctgIdx * 2 + (ctgName.second ? 0 : 1);
        if (inDegrees[i] > 0 || results[i].empty()) { continue; }

        std::stringstream nameStream;

        nameStream << prefix << nameCnt++;

        std::set<std::pair<std::size_t, bool>> connected;

        std::size_t maxLen = 0, totalLen = 0;

        combatSeq(results, *pCtgMapper, i / 2, i % 2 == 0,
                  [&](std::size_t ctgId, bool forward, std::size_t startPos)-> bool {
                      //nameStream << prefix << ctgId << "_" << (forward ? 0 : 1) << "_" << (*pContigDB)[ctgId].name() << "_" << (*pContigDB)[ctgId].size() << ";";
                      connected.emplace(ctgId, forward);
                      maxLen = std::max(maxLen, (*pContigDB)[ctgId].size());
                      totalLen += algo.seqSize(results[ctgId * 2 + (forward ? 0 : 1)]);
                      return true;
                  });

        bool isConnected = connected.size() > 1 && totalLen > maxLen * 1.05;
        bool isExtended = connected.size() == 1 && algo.seqSize(results[i]) > (*pContigDB)[ctgIdx].size() * 1.2;

        if (isConnected || isExtended) {

            const std::size_t lineSize = 70;

            std::string outPath(outDir + "/" + prefix
                                + std::to_string(i / 2) + "_" + std::to_string(i % 2) + ".fasta");

            std::string helpPath(outDir + "/" + prefix
                                + std::to_string(i / 2) + "_" + std::to_string(i % 2) + ".help");

            std::ofstream helpOut(helpPath);

            helpOut << totalLen << std::endl;
            helpOut << maxLen << std::endl;

            helpOut.close();

            std::ofstream fastaOf(outPath);

            fastaOf << ">" << nameStream.str() << std::endl;

            std::size_t cnt = 0;

            combatSeq(results, *pCtgMapper, i / 2, i % 2 == 0,
                      [&](std::size_t ctgId, bool forward, std::size_t startPos)-> bool {
                          std::cout << i << "=" << ctgId << std::endl;
                          for (auto ch : algo.seqToString(results[ctgId * 2 + (forward ? 0 : 1)], deviation,
                                                          errorRate)) {
                              fastaOf << ch;
                              if (++cnt % lineSize == 0) {
                                  fastaOf << std::endl;
                                  cnt = 0;
                              }
                          }
                          return true;
                      });

            if (cnt > 0) {
                fastaOf << std::endl;
            }
            std::cout << "Out file: " << outPath << std::endl;

            for (auto &success : connected) {
                successCtgSet.emplace((*pContigDB)[success.first].name(), success.second);
            }
        } else {
            std::cout << "Ignore output" << std::endl;
        }
    }

    return successCtgSet;
}
