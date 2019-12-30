//
// Created by godot on 19-2-18.
//

#ifndef PAGRAPH_POSITIONPROCESSER_HPP
#define PAGRAPH_POSITIONPROCESSER_HPP

#include "graph/PABruijnGraph.hpp"
#include "align/Aligner.hpp"
#include "seq/SeqInf.hpp"
#include "align/MecatAlignDatabase.hpp"
#include "align/AlignInf.hpp"
#include "PositionMapper.hpp"
#include <vector>
#include <mutex>
#include <atomic>
#include <unordered_set>
#include <memory>

class PositionProcessor {
private:
    const int _progressWidth = 70;
    const std::size_t _updateGap = 1000;

    std::shared_ptr<PABruijnGraph> _pPaGraph;

    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReadDB;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pContigDB;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReferenceDB;

    std::shared_ptr<const IAlignDatabase<AlignInf>> _pReadToContigDB;
    std::shared_ptr<const IAlignDatabase<AlignInf>> _pReadToRefDB;
    std::shared_ptr<const IAlignDatabase<AlignInf>> _pContigToRefDB;

    /* Start position */
    PositionMapper _contigPosMapper;
    PositionMapper _refPosMapper;

    bool _isFilterRef;
    std::set<std::string> _refNameWhiteList;
    std::set<std::size_t> _refIdWhiteList;

    Aligner _aligner;

    void transformPosition(
            const std::vector<std::vector<std::pair<Aligner::ref_pos_t, Aligner::ref_pos_t>>> &p1,
            std::vector<std::vector<PABruijnGraph::DualPos>> &p2);

public:
    int readToContigTopK = 3;
    int contigToRefTopK = -1;
    int readToRefTopK = 3;
    int outerSample = 1;
    int innerSample = 1;
    std::size_t posError = 10;
    double readToCtgRatio = 0.5;
    double readToRefRatio = 0.4;
    double ctgToRefRatio = 0.1;
    double ctgToRefTotalRatio = 0.3;
    std::size_t ctgToRefMinLen = 2500;
    unsigned threadNum = 8;

    explicit PositionProcessor(
            std::shared_ptr<PABruijnGraph> _pGraph,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReferenceDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToContigDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToRefDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pContigToRefDB);

    void preProcess();

    void process();

    char rv(char c);

    void setReadToCtgTopK(int topK);

    void setCtgToRefTopK(int topK);

    void setReadToRefTopK(int topK);

    void setReadToCtgRatio(double ratio);

    void setReadToRefRatio(double ratio);

    void setCtgToRefRatio(double ratio);

    void setCtgToRefTotalRatio(double ratio);

    void setOuterSample(int outerSample);

    void setInnerSample(int innerSample);

    void setPositionError(std::size_t error);

    void setThreadNum(unsigned threadNum);

    void setIsFilterRef(bool filter);

    void setRefFilter(const std::string &refName, bool accepted);

    void clearRefFilter(bool accepted);

    void setCtgFilter(const std::string &ctgName, bool forward, bool accepted);

    void clearCtgFilter(bool accepted);

    void setCtgToRefMinLen(std::size_t len);

    Aligner &getAligner();
};


#endif //PAGRAPH_POSITIONPROCESSER_HPP
