//
// Created by godot on 19-2-12.
//

#include "PAlgorithm.hpp"
#include <algorithm>

PAlgorithm::PAlgorithm(std::shared_ptr<const PABruijnGraph> g, std::shared_ptr<const ISeqDatabase<SeqInf>> pReads,
                       std::shared_ptr<const ISeqDatabase<SeqInf>> pContigs,
                       std::shared_ptr<const ISeqDatabase<SeqInf>> pReferences,
                       std::shared_ptr<const PositionMapper> pCtgMapper,
                       std::shared_ptr<const PositionMapper> pRefMapper, unsigned int threadNum) :
        _pPagraph(std::move(g)),
        _pReads(std::move(pReads)),
        _pContigs(std::move(pContigs)),
        _pReferences(std::move(pReferences)),
        _pCtgMapper(std::move(pCtgMapper)),
        _pRefMapper(std::move(pRefMapper)),
        _threadNum(threadNum),
        _globalUniqueTable(),
        _travelUniqueTable(),
        _uniqueTable(),
        _ctgGlobalPosTable(),
        _ctgTravelPosTable(),
        _ctgPosTable() {}

void PAlgorithm::filterSequence(PAlgorithm::TravelSequence &seq) {
    std::size_t windowSize = 10;

    if (seq.size() < windowSize) {
        return;
    }
    std::size_t startIdx = seq.size() - seq.size() / 90;
    std::cout << startIdx << std::endl;

    for (auto i = startIdx; i < seq.size() - windowSize + 1; ++i) {
        auto firstPos = seq[i].first.getPosition().first;
        auto secondPos = seq[std::min(seq.size(), i + windowSize) - 1].first.getPosition().first;
        if (secondPos != 0 && firstPos != 0 && secondPos < firstPos) {
            seq.resize(i + 1);
            return;
        }
    }
}

std::size_t PAlgorithm::editDistance(const std::string &str1, const std::string &str2) {
    std::vector<std::vector<std::size_t>> dp(2, std::vector<std::size_t>(str2.size() + 1, 0));
    std::size_t flag = 0U;

    for (decltype(str1.size()) j = 0; j <= str2.size(); ++j) {
        dp[flag][j] = j;
    }

    flag ^= 1U;

    for (decltype(str1.size()) i = 1; i <= str1.size(); ++i) {
        for (decltype(str2.size()) j = 0; j <= str2.size(); ++j) {
            if (j == 0) {
                dp[flag][j] = i;
            } else {
                dp[flag][j] = std::min(dp[flag ^ 1U][j] + 1, dp[flag][j - 1] + 1);
                dp[flag][j] = std::min(dp[flag][j], dp[flag ^ 1U][j - 1] + (str1[i - 1] == str2[j - 1] ? 0 : 1));
            }
        }
        flag ^= 1U;
    }

    return dp[flag ^ 1U][str2.size()];
}

bool PAlgorithm::existCtgPos(const std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table,
                             PABruijnGraph::PosType pos) {
    return pos >= table.first && pos <= table.second;
}

void PAlgorithm::resetCtgPosTable(std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table) {
    table.first = std::numeric_limits<PABruijnGraph::PosType>::max();
    table.second = 0;
}

void PAlgorithm::insertCtgPosHelper(std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table,
                                    PABruijnGraph::PosType pos) {
    if (pos == 0) {
        return;
    }

    table.first = std::min(table.first, pos);
    table.second = std::max(table.second, pos);
}

void PAlgorithm::insertCtgPosHelper(std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table,
                                    const PAlgorithm::TravelSequence &seq) {
    for (auto &it : seq) {
        auto nowCtgPos = it.first.getPosition().first;

        insertCtgPosHelper(table, nowCtgPos);
    }
}

void PAlgorithm::filterPANodes(const std::set<PABruijnGraph::PANode::UniqueType> &table,
                               std::vector<PABruijnGraph::PANode> &successors) {

    filterHelper(
            successors,
            [](const PABruijnGraph::PANode &i) -> const PABruijnGraph::PANode &
            { return i; },
            [&](const PABruijnGraph::PANode &i) { return table.count(i.getUniqueHelper()) == 0; });
}

std::int64_t
PAlgorithm::appendSeq(PAlgorithm::TravelSequence &base, const PAlgorithm::TravelSequence &tail, std::size_t deviation) {
    if (tail.empty()) {
        return 0;
    }

    std::int64_t dLen = 0;

    auto &head = tail.front();

    decltype(base.back().second) dist = _pPagraph->getKmerSize();

    while (!base.empty()
           && (base.back().first.getPosition().first == 0 ||
               head.first.getPosition().first <= base.back().first.getPosition().first)) {
        dLen -= base.back().second;
        base.pop_back();
    }

    if (!base.empty()) {
        dist = head.first.getPosition().first - base.back().first.getPosition().first;
    }

    for (auto &node : tail) {
        dLen += node.second;
        base.push_back(node);
    }

    dLen -= base[base.size() - tail.size()].second - dist;
    base[base.size() - tail.size()].second = dist;

    return dLen;
}

PAlgorithm::TravelSequence
PAlgorithm::travelSequence(std::size_t ctgIdx, bool forward, std::size_t deviation, double errorRate, double startSplit) {
    const std::size_t paNodeTopK = std::min(_threadNum, 8U);

    std::set<PABruijnGraph::PANode::UniqueType> globalUniqueTable;
    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> ctgGlobalPosTable;
    resetCtgPosTable(ctgGlobalPosTable);

    std::int64_t chosenOne = forward ? ctgIdx + 1 : -ctgIdx - 1;
    auto ctgStr = (*_pContigs)[ctgIdx].toString(forward);
    auto aNodes = _pPagraph->findAll(ctgStr);
    auto ctgLen = (*_pContigs)[ctgIdx].size();

    auto splitLen = static_cast<std::size_t>(ctgLen * startSplit);
    auto splitMin = (1 - startSplit);

    auto ctgLeft = static_cast<PABruijnGraph::PosType>(_pCtgMapper->dualToSingle(chosenOne, 0));
    auto ctgRight = static_cast<PABruijnGraph::PosType>(_pCtgMapper->dualToSingle(chosenOne, (*_pContigs)[ctgIdx].size()));
    auto revCtgLeft = static_cast<PABruijnGraph::PosType>(_pCtgMapper->dualToSingle(-chosenOne, 0));
    auto revCtgRight = static_cast<PABruijnGraph::PosType>(_pCtgMapper->dualToSingle(-chosenOne, (*_pContigs)[ctgIdx].size()));

    auto filter = [&](const PABruijnGraph::PANode &paNode, const std::pair<PABruijnGraph::PANode, int> &successor) -> bool {
        return (globalUniqueTable.count(successor.first.getUniqueHelper()) == 0)
               &&
               (successor.first.getPosition().first == 0 ||
                PABruijnGraph::isEdgeSimilar(
                        paNode.getPosition(),
                        successor.first.getPosition(),
                        successor.second,
                        deviation,
                        errorRate).first ||
                !existCtgPos(ctgGlobalPosTable, successor.first.getPosition().first))
               &&
               (successor.first.getPosition().first == 0 ||
                (successor.first.getPosition().first < revCtgLeft ||
                 successor.first.getPosition().first >= revCtgRight));
    };

    std::vector<PABruijnGraph::PANode> paNodes;

    searchPANode(
            aNodes, paNodes, true,
            [&](std::size_t originCtgPos, std::int64_t curCtgIdx, std::uint64_t curCtgPos, std::int64_t,
                std::int64_t) -> bool {

                return curCtgIdx == chosenOne &&
                       std::max(curCtgPos, originCtgPos) - std::min(curCtgPos, originCtgPos) <= deviation;
            });

    paNodes.resize(std::min(paNodes.size(), paNodeTopK));

    std::int64_t ctgStart = 0;
    for (auto &node : paNodes) {
        auto ctgPos = _pCtgMapper->singleToDual(node.getPosition().first).second;
        ctgStart = std::max(ctgStart, ctgPos);
    }

    TravelSequence travelSeq;

    TravelSequence longestSeq;
    TravelSequence leapSeq;

    std::int64_t varLen = 0;
    std::size_t countLen = 0;

    std::deque<PABruijnGraph::PosType> ctgPosQue;
    std::deque<PABruijnGraph::PosType> refPosQue;
    const std::size_t maxQueSize = 3;
    bool finalLeap = false;

    std::cout << "Ctgstart=" << ctgStart << std::endl;

    while (!paNodes.empty()) {
        //std::cout << "STEP 1" << std::endl;
        longestSeq.clear();

        std::size_t maxLen = 0;
        std::size_t chooseCtgPos = 0;
        std::size_t chooseRefPos = 0;

        bool leap = false;

        std::vector<TravelSequence> seqResult(paNodes.size());

        //std::cout << "Travel2 Start" << std::endl;

        MultiThreadTools::multiTraversal
                (paNodes, _threadNum, [&](decltype(_threadNum) t, std::size_t ii, const PABruijnGraph::PANode &paNode) {
                    //std::cout << t << " " << _pPagraph->toString(paNode) << " Start" << std::endl;
                    seqResult[ii] = graphTravel(paNode, deviation, errorRate, {ctgLeft, ctgRight}, varLen, splitLen, splitMin,
                                                *_pCtgMapper, filter);
                    //std::cout << t << " " << "Done" << std::endl;
                });

        //std::cout << "Travel2 Done" << std::endl;

        for (std::size_t i = 0; i < paNodes.size(); ++i) {
            auto &paNode = paNodes[i];
            //std::cout << "\t" << _pPagraph->toString(paNode);
            auto &seq = seqResult[i];
            auto len = seqSize(seq);
            //std::cout << " " << len << std::endl;

            leap = seq.back().first.getPosition().first != 0 &&
                   _pCtgMapper->singleToDual(seq.back().first.getPosition().first).first != chosenOne;

            if (len > maxLen || leap) {
                maxLen = len;
                longestSeq = seq;
                chooseCtgPos = _pCtgMapper->singleToDual(paNode.getPosition().first).second;
                chooseRefPos = _pRefMapper->singleToDual(paNode.getPosition().second).second;
                if (leap) {
                    break;
                }
            }
        }

        std::cout << "choose " << chooseCtgPos << ", " << chooseRefPos << std::endl;

        varLen += appendSeq(travelSeq, longestSeq, deviation);
        countLen += maxLen;

        std::cout << "\tto " << _pCtgMapper->singleToDual(longestSeq.back().first.getPosition().first).second
                  << ", " << _pRefMapper->singleToDual(longestSeq.back().first.getPosition().second).second << std::endl;

        std::cout << "\tle " << varLen << std::endl;

        //std::cout << "\t" << maxLen << std::endl;
        //std::cout << "\t" << varLen << std::endl;

        if (chooseCtgPos != 0) {
            ctgPosQue.push_back(chooseCtgPos);
            while (ctgPosQue.size() > maxQueSize) {
                ctgPosQue.pop_front();
            }
        }

        if (chooseRefPos != 0) {
            refPosQue.push_back(chooseRefPos);
            while (refPosQue.size() > maxQueSize) {
                refPosQue.pop_front();
            }
        }


        for (auto &node : longestSeq) {
            globalUniqueTable.insert(node.first.getUniqueHelper());
        }

        insertCtgPosHelper(ctgGlobalPosTable, longestSeq);

        //setAllUsed(longestSeq, true);

        bool ctgRepeatCheck = false;
        bool refRepeatCheck = false;

        if (ctgPosQue.size() >= maxQueSize) {
            auto minMax = std::minmax_element(ctgPosQue.begin(), ctgPosQue.end());
            ctgRepeatCheck = *minMax.second - *minMax.first <= 2 * deviation;
        }

        if (refPosQue.size() >= maxQueSize) {
            auto minMax = std::minmax_element(refPosQue.begin(), refPosQue.end());
            refRepeatCheck = *minMax.second - *minMax.first <= 2 * deviation;
        }

        if (ctgRepeatCheck) {
            std::cout << "REPEAT I" << std::endl;
        }

        if (refRepeatCheck) {
            std::cout << "REPEAT II" << std::endl;
        }

        if (ctgRepeatCheck || refRepeatCheck || leap) {
            if (leap) {
                std::cout << "leap" << std::endl;
                finalLeap = true;
            }
            break;
        }

        std::size_t lastCtgPos = 0;
        std::string lastCtgKmer;
        std::size_t lastRefIdx = 0;
        std::size_t lastRefPos = 0;
        std::string lastRefKmer;
        std::size_t lastLen = 0;
        bool flag1 = false, flag2 = false;

        for (auto it = travelSeq.rbegin(); (!flag1 || !flag2) && it != travelSeq.rend(); ++it) {
            auto &paNode = it->first;

            if (!flag1 && paNode.getPosition().first != 0) {
                auto dualPos = _pCtgMapper->singleToDual(paNode.getPosition().first);
                if (dualPos.first == chosenOne && dualPos.second >= 0) {
                    lastCtgPos = static_cast<std::size_t>(dualPos.second);
                    lastCtgKmer = _pPagraph->toString(paNode.getABruijnNode());
                    flag1 = true;
                }
            }

            if (!flag2) {
                lastLen += it->second;
                if (paNode.getPosition().second != 0) {
                    auto dualPos = _pRefMapper->singleToDual(paNode.getPosition().second);
                    lastRefIdx = static_cast<std::size_t>(dualPos.first);
                    lastRefPos = static_cast<std::size_t>(dualPos.second);
                    lastRefKmer = _pPagraph->toString(paNode.getABruijnNode());
                    flag2 = true;
                }
            }
        }

        paNodes.clear();

        std::string parentKmer;

        searchPANode2(
                aNodes, paNodes, false, lastCtgPos, deviation,
                [&](std::size_t, std::int64_t curCtgIdx, std::uint64_t curCtgPos, std::int64_t,
                    std::int64_t) -> bool {

                    return curCtgIdx == chosenOne
                           && std::max(curCtgPos, lastCtgPos) - std::min(curCtgPos, lastCtgPos) <= deviation;
                });
        filterPANodes(globalUniqueTable, paNodes);
        parentKmer = lastCtgKmer;


        /*if (lastRefIdx < _pReferences->size()) {
            auto refStr = (*_pReferences)[lastRefIdx].toString(true);
            auto refANodes = _pPagraph->findAll(refStr);

            searchPANode(
                    refANodes, paNodes, false,
                    [&](std::size_t cur, std::int64_t curCtgIdx, std::int64_t curCtgPos, std::int64_t curRefIdx,
                        std::int64_t curRefPos) -> bool {

                        return
                                (curCtgIdx == 0 || curCtgIdx == chosenOne) &&
                               curRefIdx == lastRefIdx &&
                               std::abs(static_cast<std::int64_t>(curRefPos - lastRefPos)) <= deviation;
                    });
            filterPANodes(globalUniqueTable, paNodes);
            parentKmer = lastRefKmer;
            //std::cout << "\t" << paNodes.size() << std::endl;
        }*/


        std::sort(paNodes.begin(), paNodes.end(),
                  [&](const PABruijnGraph::PANode &lhs, const PABruijnGraph::PANode &rhs) -> bool {
                      return editDistance(parentKmer, _pPagraph->toString((lhs.getABruijnNode()))) <
                             editDistance(parentKmer, _pPagraph->toString((rhs.getABruijnNode())));
                  });

        paNodes.resize(std::min(paNodes.size(), paNodeTopK));
    }

    if (!finalLeap) {
        filterSequence(travelSeq);
    }

    if (finalLeap) {
        auto dualPos = _pCtgMapper->singleToDual(travelSeq.back().first.getPosition().first);

        std::cout << dualPos.first << " " << dualPos.second << std::endl;

        if (static_cast<decltype(ctgIdx)>(std::abs(dualPos.first)) == ctgIdx + 1 ||
            dualPos.second >= (*_pContigs)[std::abs(dualPos.first) - 1].size() * (1 - startSplit)) {
            travelSeq.pop_back();
            std::cout << "Pump it" << std::endl;
        }
    }

    return travelSeq;
}

std::string PAlgorithm::seqToString(const PAlgorithm::TravelSequence &seq, std::size_t deviation, double errorRate,
                                    size_t ctgStartPos) {
    if (seq.empty()) {
        return "";
    }

    std::string str;

    str.append(_pPagraph->toString(seq[0].first.getABruijnNode()));

    auto kmerSize = static_cast<int>(_pPagraph->getKmerSize());

    auto firstPos = static_cast<std::size_t>(_pCtgMapper->singleToDual(seq.front().first.getPosition().first).second);

    for (std::size_t i = 1; i < seq.size(); ++i) {
        auto &prevItem = seq[i - 1];
        auto &nowItem = seq[i];

        auto similar = PABruijnGraph::isEdgeSimilar(prevItem.first.getPosition(), nowItem.first.getPosition(), nowItem.second, deviation, errorRate);

        auto useCtg = similar.first;

        if (!similar.first && !similar.second) {
            auto posSimilar = PABruijnGraph::isPosSimilar(prevItem.first.getPosition(), nowItem.first.getPosition(), deviation);
            useCtg = posSimilar.first;
        }

        auto &selectedRefDB = useCtg ? *_pContigs : *_pReferences;
        auto &selectedMapper = useCtg ? *_pCtgMapper : *_pRefMapper;
        auto selectedPrevPos = useCtg ? prevItem.first.getPosition().first : prevItem.first.getPosition().second;
        auto selectedNowPos = useCtg ? nowItem.first.getPosition().first : nowItem.first.getPosition().second;

        auto startPos = selectedMapper.singleToDual(selectedPrevPos);
        auto endPos = selectedMapper.singleToDual(selectedNowPos);

        //auto kmerDist = useCtg && endPos.second > startPos.second ?
        //        static_cast<decltype(nowItem.second)>(endPos.second - startPos.second) : nowItem.second;
        auto kmerDist = nowItem.second;
        auto posDist = endPos.second - startPos.second;

        auto selectedRefIdx = std::abs(endPos.first) - 1;
        auto selectedRefForward = endPos.first > 0;

        double moveLen = posDist * 1.0 / kmerDist;
        double refNow = startPos.second + kmerSize;

        auto kmer = _pPagraph->toString(nowItem.first.getABruijnNode());

        for (int j = 0; j < kmerDist; ++j) {
            if (kmerSize - kmerDist + j >= 0) {
                str.push_back(kmer[kmerSize - kmerDist + j]);
            } else {
                auto refPos = static_cast<std::size_t>(std::round(refNow));
                str.push_back(static_cast<char>(std::tolower(selectedRefDB[selectedRefIdx].quickBaseAt(refPos, selectedRefForward))));
            }

            refNow += moveLen;
        }
    }

    return str.substr(ctgStartPos >= firstPos ? ctgStartPos - firstPos : 0);
}

std::size_t PAlgorithm::seqSize(const PAlgorithm::TravelSequence &seq) {
    std::size_t size = 0;
    for (auto &node : seq) {
        size += node.second;
    }
    return size;
}
