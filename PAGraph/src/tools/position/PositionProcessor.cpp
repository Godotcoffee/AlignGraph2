//
// Created by godot on 19-2-18.
//

#include <chrono>

#include "PositionProcessor.hpp"

PositionProcessor::PositionProcessor(std::shared_ptr<PABruijnGraph> _pGraph,
                                     std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
                                     std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
                                     std::shared_ptr<const ISeqDatabase<SeqInf>> pReferenceDB,
                                     std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToContigDB,
                                     std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToRefDB,
                                     std::shared_ptr<const IAlignDatabase<AlignInf>> pContigToRefDB)
        : _pPaGraph(std::move(_pGraph)),
          _pReadDB(std::move(pReadDB)),
          _pContigDB(std::move(pContigDB)),
          _pReferenceDB(std::move(pReferenceDB)),
          _pReadToContigDB(std::move(pReadToContigDB)),
          _pReadToRefDB(std::move(pReadToRefDB)),
          _pContigToRefDB(std::move(pContigToRefDB)),
          _contigPosMapper(*_pContigDB),
          _refPosMapper(*_pReferenceDB),
          _isFilterRef(false),
          _refNameWhiteList(),
          _refIdWhiteList(),
          _aligner(
                  _pReadDB,
                  _pContigDB,
                  _pReferenceDB,
                  _pReadToContigDB,
                  _pReadToRefDB,
                  _pContigToRefDB
          ) {}

void PositionProcessor::transformPosition(
        const std::vector<std::vector<std::pair<Aligner::ref_pos_t, Aligner::ref_pos_t>>> &p1,
        std::vector<std::vector<PABruijnGraph::DualPos>> &p2) {

    p2.clear();
    p2.resize(p1.size());

    for (std::size_t i = 0; i < p1.size(); ++i) {
        for (auto &pos : p1[i]) {
            //if (!_isFilterRef || (pos.second.first == 0 || _refIdWhiteList.count(pos.second.first - 1) > 0)) {
                p2[i].emplace_back(
                        static_cast<PABruijnGraph::PosType>(_contigPosMapper.dualToSingle(pos.first.first,
                                                                                          pos.first.second)),
                        static_cast<PABruijnGraph::PosType>(_refPosMapper.dualToSingle(pos.second.first,
                                                                                       pos.second.second)));
            //}
        }
    }
}

void
PositionProcessor::preProcess() {
    auto &aligner = _aligner;

    aligner.simpleAlign();

    /*aligner.clearCtgFilter(false);
    for (auto &cc : ctgCov) {
        auto ctgLen = (*_pContigDB)[_pContigDB->seqId(cc.first)].size();
        auto cov = std::max(
                1.0 * cc.second.first / ctgLen,
                1.0 * cc.second.second / ctgLen
        );
        if (cov >= ctgToRefTotalRatio) {
            aligner.setCtgFilter(cc.first, true);
        }
    }*/

    aligner.addExtraPosition();
}

void
PositionProcessor::process() {
    auto &aligner = _aligner;
    auto &readDB = *_pReadDB;
    auto &paGraph = *_pPaGraph;

    auto totalRead = _pReadDB->size();

    std::atomic_size_t readCnt;
    std::mutex lock;

    auto functor = [&](
            std::size_t idx,
            std::vector<std::vector<std::pair<Aligner::ref_pos_t, Aligner::ref_pos_t>>> &fp,
            std::vector<std::vector<std::pair<Aligner::ref_pos_t, Aligner::ref_pos_t>>> &bp,
            bool fwdUseful,
            bool revUseful,
            unsigned threadId) {

        if (fwdUseful) {
            std::vector<std::vector<PABruijnGraph::DualPos>> posAdapter;
            auto forward = readDB[idx].toString(true);
            transformPosition(fp, posAdapter);
            paGraph.addPositionAndEdge(forward, posAdapter, outerSample);
        }

        if (revUseful) {
            std::vector<std::vector<PABruijnGraph::DualPos>> posAdapter;
            auto revc = readDB[idx].toString(false);
            transformPosition(bp, posAdapter);
            paGraph.addPositionAndEdge(revc, posAdapter, outerSample);
        }

        if (++readCnt % _updateGap == 0) {
            std::lock_guard<std::mutex> lockGuard(lock);
            MyTools::progress(readCnt * 1.0 / totalRead, _progressWidth, std::cerr);
        }
    };

    std::cout << "[PositionProcessor] Running read to contig..." << std::endl;
    MyTools::progress(0.0, _progressWidth, std::cerr);
    readCnt = 0;

    aligner.runReadToCtg(functor, threadNum);

    MyTools::progress(1.0, _progressWidth, std::cerr);
    std::cout << std::endl;

    std::cout << "\tmerge edge = "  << paGraph.mergeEdge(false, threadNum)                       << std::endl;
    std::cout << "\ttotal pos = "   << paGraph.totalPosition(threadNum)                          << std::endl;
    std::cout << "\tmerge pos = "   << paGraph.mergeKmerPosition(posError, false, threadNum)     << std::endl;

    std::cout << std::endl;

    MyTools::progress(0.0, _progressWidth, std::cerr);
    readCnt = 0;

    aligner.runReadToRef(functor, threadNum);

    MyTools::progress(1.0, _progressWidth, std::cerr);
    std::cout << std::endl;

    std::cout << "\tmerge edge = "  << paGraph.mergeEdge(true, threadNum)                       << std::endl;
    std::cout << "\ttotal pos = "   << paGraph.totalPosition(threadNum)                         << std::endl;
    std::cout << "\tmerge pos = "   << paGraph.mergeKmerPosition(posError, true, threadNum)     << std::endl;

    paGraph.sortKmerPosition(threadNum);

    paGraph.resetUsedFlag(threadNum);

    std::cout << "[PositionProcessor] Done!" << std::endl;

    //return ctgCov;
}

char PositionProcessor::rv(char c) {
    switch(c) {
        case 'A': case 'a': default:
            return 'T';
        case 'C': case 'c':
            return 'G';
        case 'G': case 'g':
            return 'C';
        case 'T': case 't':
            return 'A';
    }
}

void PositionProcessor::setReadToCtgTopK(int topK) {
    this->readToContigTopK = topK;
    _aligner.setReadToCtgTopK(this->readToContigTopK);
}

void PositionProcessor::setCtgToRefTopK(int topK) {
    this->contigToRefTopK = topK;
    _aligner.setCtgToRefTopK(this->contigToRefTopK);
}

void PositionProcessor::setReadToRefTopK(int topK) {
    this->readToRefTopK = topK;
    _aligner.setReadToRefTopK(this->readToRefTopK);

}

void PositionProcessor::setReadToCtgRatio(double ratio) {
    this->readToCtgRatio = ratio;
    _aligner.setReadToCtgRatio(this->readToCtgRatio);
}

void PositionProcessor::setReadToRefRatio(double ratio) {
    this->readToRefRatio = ratio;
    _aligner.setReadToRefRatio(this->readToRefRatio);
}

void PositionProcessor::setCtgToRefRatio(double ratio) {
    this->ctgToRefRatio = ratio;
    _aligner.setCtgToRefRatio(this->ctgToRefRatio);
}

void PositionProcessor::setCtgToRefTotalRatio(double ratio) {
    this->ctgToRefTotalRatio = ratio;
    _aligner.setCtgToRefTotalRatio(this->ctgToRefTotalRatio);
}

void PositionProcessor::setOuterSample(int outerSample) {
    this->outerSample = outerSample;
}

void PositionProcessor::setInnerSample(int innerSample) {
    this->innerSample = innerSample;
}

void PositionProcessor::setPositionError(std::size_t error) {
    this->posError = error;
}

void PositionProcessor::setCovFilter(std::size_t cov) {
    _aligner.setCovFilter(cov);
}

void PositionProcessor::setThreadNum(unsigned threadNum) {
    this->threadNum = threadNum;
}

void PositionProcessor::setIsFilterRef(bool filter) {
    this->_isFilterRef = filter;
}

void PositionProcessor::setRefFilter(const std::string &refName, bool accepted) {
    _aligner.setRefFilter(refName, accepted);
}

void PositionProcessor::clearRefFilter(bool accepted) {
    _aligner.clearRefFilter(accepted);
}

void PositionProcessor::setCtgFilter(const std::string &refName, bool forward, bool accepted) {
    _aligner.setCtgFilter(refName, forward, accepted);
}

void PositionProcessor::clearCtgFilter(bool accepted) {
    _aligner.clearCtgFilter(accepted);
}

void PositionProcessor::setCtgToRefMinLen(std::size_t len) {
    this->ctgToRefMinLen = len;
    _aligner.setCtgToRefMinLen(this->ctgToRefMinLen);
}

Aligner &PositionProcessor::getAligner() {
    return _aligner;
}
