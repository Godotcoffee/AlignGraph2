//
// Created by alice on 18-10-3.
//

#include "PABruijnGraph.hpp"
#include <atomic>
#include "kmer/KmerHelper.hpp"
#include "tools/MyTools.hpp"

PABruijnGraph::PABruijnGraph(const IKmerIterator &kmerIt, unsigned threadNum) :
        _kmerSize(kmerIt.kSize()),
        _indexMask((1UL << (_kmerSize * 2)) - 1) {

    std::vector<KmerIndex> idxAbundance;

    if (threadNum == 1) {
        kmerIt.iterate([&](std::uint64_t kmer) {
            idxAbundance.emplace_back(kmer);
        }, threadNum);
    } else {
        std::mutex mutex;
        kmerIt.iterate([&](const std::uint64_t kmer) {
            {
                std::lock_guard<std::mutex> lockGuard(mutex);
                idxAbundance.emplace_back(kmer);
            }
        }, threadNum);
    }

    std::sort(idxAbundance.begin(), idxAbundance.end());
    idxAbundance.erase(std::unique(idxAbundance.begin(), idxAbundance.end()), idxAbundance.end());

    for (auto &p : idxAbundance) {
        _kmerIndexArr.push_back(p);
    }

    for (std::size_t i = 0; i < _kmerIndexArr.size(); ++i) {
        _kmerIndexMap[_kmerIndexArr[i]] = i;
    }

    //_pDenseHashTable = std::make_unique<DenseHashTable>(_kmerIndexArr.size());
    std::cout << "Hash Table Size = " << _kmerIndexArr.size() << std::endl;
    std::cout << "Hash Table Memory = " << _kmerIndexArr.size() * sizeof(KMerNode) / 1024 << "KB" << std::endl;
    _pDenseHashTable = MyTools::make_unique<DenseHashTable>(_kmerIndexArr.size());
}

char PABruijnGraph::revAcgt(char ch) const {
    switch (ch) {
        case 'A': case 'a':
            return 'T';
        case 'C': case 'c':
            return 'G';
        case 'G': case 'g':
            return 'C';
        case 'T': case 't': default:
            return 'A';
    }
}

unsigned PABruijnGraph::acgt(char ch) const {
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

PABruijnGraph::KmerIndex PABruijnGraph::calcKmerCode(const std::string &kmer, int start) const {
    std::vector<KmerIndex> codes;

    KmerHelper::kmer2Code(codes, kmer, _kmerSize);
    return !codes.empty() ? codes.front() : 0;
}

std::string PABruijnGraph::calcKmerString(PABruijnGraph::KmerIndex kmerCode) const {
    return KmerHelper::code2Kmer(kmerCode, _kmerSize);
}

std::pair<PABruijnGraph::ANode, bool> PABruijnGraph::fromCode(PABruijnGraph::KmerIndex code) const {
    auto check = searchDenseIndex(code);
    if (!check.second) {
        return {ANode(), false};
    }
    return {ANode((*_pDenseHashTable)[check.first], check.first), true};
}

std::pair<PABruijnGraph::KmerIndex, bool> PABruijnGraph::searchDenseIndex(PABruijnGraph::KmerIndex code) const {
    //auto ptr = std::lower_bound(_kmerIndexArr.begin(), _kmerIndexArr.end(), idx);
    //return std::make_pair(ptr - _kmerIndexArr.begin(), ptr != _kmerIndexArr.end() && *ptr == idx);
    auto it = _kmerIndexMap.find(code);
    auto exist = it != _kmerIndexMap.end();
    return {exist ? it->second : -1, exist};
}

std::vector<PABruijnGraph::KmerIndex> PABruijnGraph::buildCode(const std::string &seq) const {
    //KmerIndex rcCode = 0;
    //std::vector<KmerDualIndex> codes;
    std::vector<KmerIndex> codes;

    KmerHelper::kmer2Code(codes, seq, _kmerSize);

    /*auto halfSize = codes.size() / 2;
    for (decltype(halfSize) i = 0; i < halfSize; ++i) {
        std::swap(codes[i].second, codes[codes.size() - 1 - i].second);
    }*/

    // 不确定编译器是否会进行RVO
    return codes;
}

bool PABruijnGraph::preCheckPos(const PABruijnGraph::DualPos &pos1, const PABruijnGraph::DualPos &pos2) const {
    return !(pos1.first == 0 && pos1.second == 0) &&
           !(pos2.first == 0 && pos2.second == 0) &&
           !(pos1.first == 0 && pos2.second == 0) &&
           !(pos1.second == 0 && pos2.first == 0);

    /*if (pos1.first == 0 && pos1.second == 0) {
        return false;
    }
    if (pos2.first == 0 && pos2.second == 0) {
        return false;
    }
    if (pos1.first == 0 && pos2.second == 0) {
        return false;
    }
    if (pos1.second == 0 && pos2.first == 0) {
        return false;
    }
    return true;*/
}

PABruijnGraph::MatchGrade
PABruijnGraph::checkPosition(const PABruijnGraph::DualPos &pos1, const PABruijnGraph::DualPos &pos2,
                             PABruijnGraph::PosType dist, PABruijnGraph::PosType deviation, double errorRate) {
    //if (!preCheckPos(pos1, pos2)) {
    //    return Oops;
    //}

    auto state = isEdgeSimilar(pos1, pos2, dist, deviation, errorRate);

    auto state1 = state.first;  // Contig
    auto state2 = state.second; // Reference

    state1 = state1 || std::abs(1.0 - ((pos2.first - pos1.first)   * 1.0 / dist)) <= errorRate;
    state2 = state2 || std::abs(1.0 - ((pos2.second - pos1.second) * 1.0 / dist)) <= errorRate;

    if (pos1.first == 0 || pos2.first == 0) {
        return state2 ? (pos2.first != 0 ? Excellent : (pos1.first != 0 ? Skip : Good)) : Oops;
    }
    if (pos1.second == 0 || pos2.second == 0) {
        return state1 ? (pos2.second != 0 ? Excellent : Good) : Oops;
    }
    return (state1 && state2) ? Amazing : (state1 ? Excellent : (state2 ? Skip: Oops));
}

void
PABruijnGraph::searchSuccessors(const PABruijnGraph::PANode &parent, std::vector<std::pair<PABruijnGraph::PANode, int>> &result,
                                int deviation, double errorRate) const {
    auto &hashTable = *_pDenseHashTable;

    auto rootPos = parent.getPosition();
    auto &parentANode = parent.getABruijnNode();

    auto &node = hashTable[parentANode.getIndex()];
    auto &tos = node.getAllChild();

    std::vector<std::pair<PANode, int>> amazing;    // best
    std::vector<std::pair<PANode, int>> excellent;  // better
    std::vector<std::pair<PANode, int>> great;      // good

    for (auto toIdx : tos) {
        ANode aNode(hashTable[toIdx.first], toIdx.first);
        for (decltype(aNode.size()) j = 0; j < aNode.size(); ++j) {
            if (aNode.isUsed(j)) { continue; }

            auto pos = aNode[j];

            auto check = checkPosition(rootPos, pos, toIdx.second, deviation, errorRate);

            if (check != Oops) {
                result.emplace_back(PANode(aNode, j), toIdx.second);
            }
        }
    }

}

PABruijnGraph::PosType PABruijnGraph::reversePosition(PABruijnGraph::PosType pos) const {
    if (pos == 0) { return 0; }
    return -pos + 1 - (PosType) _kmerSize;
}

PABruijnGraph::DualPos PABruijnGraph::reversePositions(const PABruijnGraph::DualPos &pos) const {
    DualPos newPos = pos;
    if (pos.first != 0) {
        newPos.first = reversePosition(pos.first);
    }
    if (pos.second != 0) {
        newPos.second = reversePosition(pos.second);
    }
    return newPos;
}

PABruijnGraph::DualPos PABruijnGraph::buildPosition(const PABruijnGraph::DualPos &pos, bool revComp) const {
    return revComp ? reversePositions(pos) : pos;
}

void PABruijnGraph::addPosition(PABruijnGraph::KmerIndex idx, const PABruijnGraph::DualPos &pos, bool useMutex) {
    auto &node = (*_pDenseHashTable)[idx];
    node.addPosition(pos, useMutex);
}

void PABruijnGraph::addPosition(PABruijnGraph::KmerIndex idx, const std::vector<PABruijnGraph::DualPos> &pos, bool useMutex) {
    auto &node = (*_pDenseHashTable)[idx];
    node.addPosition(pos, useMutex);
}

void PABruijnGraph::addEdge(PABruijnGraph::KmerIndex idx, PABruijnGraph::KmerIndex child, int step, bool useMutex) {
    auto &node = (*_pDenseHashTable)[idx];
    node.addChild({child, step}, useMutex);
}

bool PABruijnGraph::checkIndex(PABruijnGraph::KmerIndex idx) const {
    return idx >= 0 && idx < _pDenseHashTable->size();
}

void PABruijnGraph::addPositionAndEdge(const std::string &seq, const std::vector<std::vector<PABruijnGraph::DualPos>> &fwd,
                                       int outerSample) {
    std::vector<std::pair<KmerIndex, std::size_t>> sampleResult;
    sampleSequence(sampleResult, seq, outerSample, [&](std::size_t ii) {
        return !fwd[ii].empty();
    });

    for (auto &sample : sampleResult) {
        auto kmerIndex = sample.first;
        auto posIndex = sample.second;
        addPosition(kmerIndex, fwd[posIndex]);
    }

    for (decltype(sampleResult.size()) i = 1; i < sampleResult.size(); ++i) {
        auto &from = sampleResult[i - 1];
        auto &to = sampleResult[i];
        auto distance = (int) (to.second - from.second);
        addEdge(from.first, to.first, distance);
    }
}

std::size_t PABruijnGraph::mergeKmerPosition(std::size_t error, bool shrink, unsigned int threadNum) {
    std::atomic_size_t count(0);
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) -> void {
                count += node.mergePosition([&](const DualPos &lhs, const DualPos &rhs) -> bool {
                    auto similar = isPosSimilar(lhs, rhs, error);
                    similar.first = similar.first || (lhs.first == 0 && rhs.first == 0);
                    similar.second = similar.second || (lhs.second == 0 && rhs.second == 0);
                    return similar.first && similar.second;
                }, shrink);
            });

    return count;
}

void PABruijnGraph::sortKmerPosition(unsigned int threadNum) {
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) -> void {
                node.sortPosition();
            });
}

std::size_t PABruijnGraph::mergeEdge(bool shrink, unsigned int threadNum) {
    std::atomic_size_t count(0);
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) -> void {
                count += node.mergeChild([](const AdjEdge &lhs, const AdjEdge &rhs) -> bool {
                    return lhs == rhs;
                }, shrink);
            });

    return count;
}

void PABruijnGraph::resetUsedFlag(unsigned int threadNum) {
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) {
                node.resetUsedFlag();
            });

}

void PABruijnGraph::resetAllNodes(unsigned int threadNum) {
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) {
                node.resetAll(false);
            });
}

std::size_t PABruijnGraph::totalPosition(unsigned int threadNum) const {
    std::atomic_size_t count(0);
    auto &hashTable = *_pDenseHashTable;

    MultiThreadTools::multiTraversal
            (hashTable, threadNum, [&](decltype(threadNum), decltype(hashTable.size()), KMerNode &node) {
                count += node.getAllPositions().size();
            });

    return count;
}

std::size_t PABruijnGraph::getKmerSize() const {
    return _kmerSize;
}

std::pair<PABruijnGraph::ANode, bool> PABruijnGraph::find(const std::string &kmer) const {
    auto kmerCode = calcKmerCode(kmer);
    return fromCode(kmerCode);
}

std::vector<std::pair<PABruijnGraph::ANode, size_t>> PABruijnGraph::findAll(const std::string &seq) const {
    std::vector<std::pair<ANode, std::size_t>> aNodes;
    auto index = buildCode(seq);

    auto len = seq.size();
    if (len >= _kmerSize) {
        for (std::size_t i = 0; i < len - _kmerSize + 1; ++i) {
            auto node = fromCode(index[i]);
            if (node.second) {
                aNodes.emplace_back(node.first, i);
            }
        }
    }
    return aNodes;
}

std::string PABruijnGraph::toString(const ANode &aNode) const {
    if (checkIndex(aNode.getIndex())) {
        return calcKmerString(_kmerIndexArr[aNode.getIndex()]);
    }
    return "";
}

std::string PABruijnGraph::toString(const PANode &paNode) const {
    std::stringstream ss;
    auto pos = paNode.getPosition();

    ss << toString(paNode.getABruijnNode()) << "," << pos.first << "," << pos.second << "," << paNode.getAbundance();
    return ss.str();
}

void
PABruijnGraph::successors(std::vector<std::pair<PANode, int>> &results, const PANode &paNode, PosType error, double errorRate) const {
    searchSuccessors(paNode, results, error, errorRate);
}

std::size_t PABruijnGraph::availableKmerNumber() const {
    return _kmerIndexArr.size();
}

std::pair<bool, bool> PABruijnGraph::isPosSimilar(const DualPos &lhs, const DualPos &rhs, std::size_t deviation) {
    auto state1 = ((lhs.first  != 0 && rhs.first  != 0 && std::max(lhs.first, rhs.first) - std::min(lhs.first, rhs.first)  <= deviation));
    auto state2 = ((lhs.second != 0 && rhs.second != 0 && std::max(lhs.second, rhs.second) - std::min(lhs.second, rhs.second) <= deviation));
    return {state1, state2};
}

std::pair<bool, bool>
PABruijnGraph::isEdgeSimilar(const DualPos &lhs, const DualPos &rhs, int dist, std::size_t deviation,
                             double errorRate) {
    auto tmpPos = DualPos(lhs.first != 0 ? lhs.first + dist : 0, lhs.second != 0 ? lhs.second + dist : 0);

    auto state = isPosSimilar(tmpPos, rhs, deviation);
    auto state1 = state.first;
    auto state2 = state.second;

    state1 = state1 || (lhs.first != 0 && rhs.first != 0 &&
                        std::abs(1.0 - ((rhs.first  - lhs.first)  * 1.0 / dist)) <= errorRate);
    state2 = state2 || (lhs.second != 0 && rhs.second != 0 &&
                        std::abs(1.0 - ((rhs.second - lhs.second) * 1.0 / dist)) <= errorRate);

    return {state1, state2};
}
