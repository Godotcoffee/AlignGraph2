#include <utility>

//
// Created by godot on 19-2-12.
//

#ifndef PAGRAPH_PALGORITHM_HPP
#define PAGRAPH_PALGORITHM_HPP

#include "PABruijnGraph.hpp"
#include "align/Aligner.hpp"
#include "position/PositionMapper.hpp"
#include <vector>
#include <unordered_set>
#include <string>
#include <fstream>
#include <memory>
#include <limits>
#include <set>

class PAlgorithm {
public:
    typedef std::vector<std::pair<PABruijnGraph::PANode, int>> TravelSequence;
private:
    std::shared_ptr<const PABruijnGraph> _pPagraph;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReads;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pContigs;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReferences;

    std::shared_ptr<const PositionMapper> _pCtgMapper;
    std::shared_ptr<const PositionMapper> _pRefMapper;

    unsigned int _threadNum;

    enum NodeStatus {
        End,
        Branch,
        Limit,
        Leap
    };

    std::set<PABruijnGraph::PANode::UniqueType> _globalUniqueTable;
    std::set<PABruijnGraph::PANode::UniqueType> _travelUniqueTable;
    std::set<PABruijnGraph::PANode::UniqueType> _uniqueTable;

    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> _ctgGlobalPosTable;
    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> _ctgTravelPosTable;
    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> _ctgPosTable;

    static void filterSequence(TravelSequence &seq);

public:
    static std::size_t editDistance(const std::string &str1, const std::string &str2);
private:
    template<typename Filter>
    void classifySuccessors(
            std::vector<std::pair<PABruijnGraph::PANode, int>> &results,
            const PABruijnGraph::PANode &paNode, PABruijnGraph::PosType deviation, double errorRate,
            std::pair<std::int64_t, std::int64_t> ctgRange, bool canLeap, double leapMin,
            const PositionMapper &ctgMapper, Filter filter);

    static bool existCtgPos(const std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table, PABruijnGraph::PosType pos);

    static void resetCtgPosTable(std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table);

    static void insertCtgPosHelper(
            std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table,
            PABruijnGraph::PosType pos);

    static void insertCtgPosHelper(std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> &table,
                            const TravelSequence &seq);

    template<typename Selector, typename Filter, typename T>
    static void filterHelper(std::vector<T> &contents, Selector s, Filter f);

    template<typename Filter>
    static void filterSuccessors(const PABruijnGraph::PANode &paNode, std::vector<std::pair<PABruijnGraph::PANode, int>> &successors, Filter f);

    static void filterPANodes(
            const std::set<PABruijnGraph::PANode::UniqueType> &table,
            std::vector<PABruijnGraph::PANode> &successors);

    template<typename ParentFilter>
    NodeStatus walkStraight(
            const std::pair<PABruijnGraph::PANode, int> &startNode,
            std::vector<std::pair<PABruijnGraph::PANode, int>> &path,
            int deviation, double errorRate, std::pair<std::int64_t, std::int64_t> ctgRange,
            std::size_t hasSize, std::size_t splitSize, double splitMin, const PositionMapper &ctgMapper,
            ParentFilter parentFilter, std::size_t limitation = 0);

    template<typename ParentFilter>
    TravelSequence graphTravel(
            const PABruijnGraph::PANode &startNode, int deviation, double errorRate,
            std::pair<std::int64_t, std::int64_t> ctgRange,
            std::size_t hasSize, std::size_t splitSize, double splitMin, const PositionMapper &ctgMapper, ParentFilter parentFilter);
    template<typename Functor>
    void searchPANode(
            std::vector<std::pair<PABruijnGraph::ANode, size_t>> &aNodes,
            std::vector<PABruijnGraph::PANode> &result,
            bool onlyFirst,
            Functor f);

    template<typename Functor>
    void searchPANode2(
            std::vector<std::pair<PABruijnGraph::ANode, size_t>> &aNodes,
            std::vector<PABruijnGraph::PANode> &result,
            bool onlyFirst, std::size_t pos, std::size_t deviation,
            Functor f);

    std::int64_t appendSeq(TravelSequence &base, const TravelSequence &tail, std::size_t deviation);
public:
    explicit PAlgorithm(
            std::shared_ptr<const PABruijnGraph> g,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReads,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pContigs,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReferences,
            std::shared_ptr<const PositionMapper> pCtgMapper,
            std::shared_ptr<const PositionMapper> pRefMapper,
            unsigned int threadNum);

    TravelSequence travelSequence(std::size_t ctgIdx, bool forward, std::size_t deviation, double errorRate,
                                  double startSplit);

    std::string seqToString(const TravelSequence &seq, std::size_t deviation,
                            double errorRate, std::size_t ctgStartPos = 0);

    static std::size_t seqSize(const TravelSequence &seq);
};

#include "PAlgorithm.tcc"

#endif //PAGRAPH_PALGORITHM_HPP
