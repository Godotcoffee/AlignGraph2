//
// Created by alice on 18-10-3.
//

#ifndef PAGRAPH_ABRUIJNGRAPH_HPP
#define PAGRAPH_ABRUIJNGRAPH_HPP

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <thread>
#include <algorithm>
#include <memory>
#include <tools/MyTools.hpp>
#include "kmer/IKmerIterator.hpp"
#include "seq/ISeqDatabase.hpp"
#include "seq/SeqInf.hpp"
#include "node/KMerAdjNode.hpp"
#include "node/ABruijnNode.hpp"
#include "node/PABruijnNode.hpp"
#include "thread/MultiThreadTools.hpp"

class PABruijnGraph {
public:
    typedef std::uint32_t PosType;
    typedef std::pair<PosType, PosType> DualPos;
    typedef std::uint16_t CountType;

private:
    typedef std::uint64_t KmerIndex;
    //typedef std::pair<KmerIndex, KmerIndex> KmerDualIndex;
    typedef std::pair<std::uint64_t, int> AdjEdge;
    typedef KMerAdjNode<DualPos, AdjEdge> KMerNode;

public:
    typedef ABruijnNode<DualPos, AdjEdge, CountType> ANode;
    typedef PABruijnNode<DualPos, AdjEdge, CountType> PANode;

    enum MatchGrade {
        Oops, Skip, Good, Excellent, Amazing
    };
private:
    std::vector<KmerIndex> _kmerIndexArr;
    std::unordered_map<KmerIndex, KmerIndex> _kmerIndexMap;

    typedef std::vector<KMerNode> DenseHashTable;

    std::unique_ptr<DenseHashTable> _pDenseHashTable;
    //DenseHashTable *_denseHashTable;    // For late initialization. Released in destructor.

    std::size_t _kmerSize;
    std::size_t _indexMask;

    char revAcgt(char ch) const;

    unsigned acgt(char ch) const;

    KmerIndex calcKmerCode(const std::string &kmer, int start = 0) const;

    std::string calcKmerString(KmerIndex kmerCode) const;

    std::pair<ANode, bool> fromCode(KmerIndex code) const;

    std::pair<KmerIndex, bool> searchDenseIndex(KmerIndex code) const;

    std::vector<KmerIndex> buildCode(const std::string &seq) const;

    bool preCheckPos(const DualPos &pos1, const DualPos &pos2) const;

    void searchSuccessors(const PANode &parent, std::vector<std::pair<PANode, int>> &result, int deviation, double errorRate) const;

    PosType reversePosition(PosType pos) const;

    DualPos reversePositions(const DualPos &pos) const;

    DualPos buildPosition(const DualPos &pos, bool revComp) const;

    void addPosition(KmerIndex idx, const DualPos &pos, bool useMutex = true);

    void addPosition(KmerIndex idx, const std::vector<DualPos> &pos, bool useMutex = true);

    void addEdge(KmerIndex idx, KmerIndex child, int step, bool useMutex = true);

    bool checkIndex(KmerIndex idx) const;

    template<typename T>
    void sampleSequence(std::vector<std::pair<KmerIndex, std::size_t>> &result, const std::string &seq, std::size_t outerSample, T functor) const;
public:
    explicit PABruijnGraph(const IKmerIterator &kmerIt, unsigned threadNum = 8);

    //void generateEdges(const ISeqDatabase<SeqInf> &seqDB, int outerSample, bool reverse = true, unsigned threadNum = 8);

    //void addSeqPosition(const std::string &seq, const std::vector<std::vector<DualPos>> &fwd, int outerSample);

    void addPositionAndEdge(const std::string &seq, const std::vector<std::vector<PABruijnGraph::DualPos>> &fwd, int outerSample);

    std::size_t mergeKmerPosition(std::size_t error, bool shrink = false, unsigned threadNum = 8);

    void sortKmerPosition(unsigned threadNum = 8);

    std::size_t mergeEdge(bool shrink = false, unsigned threadNum = 8);

    void resetUsedFlag(unsigned threadNum = 8);

    void resetAllNodes(unsigned threadNum = 8);

    std::size_t totalPosition(unsigned threadNum = 8) const;

    std::size_t getKmerSize() const;

    std::pair<ANode, bool> find(const std::string &kmer) const;

    std::vector<std::pair<ANode, std::size_t>> findAll(const std::string &seq) const;

    std::string toString(const ANode &aNode) const;

    std::string toString(const PANode &paNode) const;

    void successors(std::vector<std::pair<PANode, int>> &results, const PANode &paNode, PosType error = 20, double errorRate = 0.2) const;

    std::size_t availableKmerNumber() const;

    static std::pair<bool, bool> isPosSimilar(const DualPos &lhs, const DualPos &rhs, std::size_t deviation = 20);

    static std::pair<bool, bool>
    isEdgeSimilar(const DualPos &lhs, const DualPos &rhs, int dist, std::size_t deviation = 20,
                  double errorRate = 0.2);

    static MatchGrade
    checkPosition(const DualPos &pos1, const DualPos &pos2, PosType dist, PosType deviation, double errorRate);
};

#include "PABruijnGraph.tcc"

#endif //PAGRAPH_ABRUIJNGRAPH_HPP
