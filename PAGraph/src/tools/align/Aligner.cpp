//
// Created by alice on 18-11-7.
//

#include "Aligner.hpp"
#include <map>

void Aligner::init(bool debug) {
    if (!_isInit) {
        if (debug) {
            std::cout << "=== ARGS ===" << std::endl;
            std::cout << "ctg to ref ratio = " << _ctgToRefRatio << std::endl;
        }

        if (debug) { std::cout << "read number: " << _pReadDB->size() << std::endl; }
        if (debug) { std::cout << "contig number: " << _pContigDB->size() << std::endl; }
        if (debug) { std::cout << "reference number: " << _pReferenceDB->size() << std::endl; }

        if (debug) { std::cout << "merge align information" << std::endl; }
        mergeAlignInf();

        if (debug) { std::cout << "simple align" << std::endl; }
        simpleAlign();

        if (debug) { std::cout << "add extra pos" << std::endl; }
        addExtraPosition();

        _isInit = true;
    }
}

void Aligner::mergeAlignInfHelper(std::vector<std::vector<Aligner::ExAlignInf>> &align,
                                  const IAlignDatabase<AlignInf> &alignDB,
                                  const ISeqDatabase<SeqInf> &queryDB, const ISeqDatabase<SeqInf> &refDB) {
    align.clear();
    align.resize(queryDB.size());

    std::size_t total = alignDB.size();

    for (decltype(total) i = 0; i < total; ++i) {
        auto &alignInf = alignDB[i];
        auto &queryName = alignInf.getQueryName();
        auto &refName = alignInf.getRefName();

        if (queryDB.containsName(queryName) && refDB.containsName(refName)) {
            auto queryIdx = queryDB.seqId(queryName);
            auto refIdx = refDB.seqId(refName);

            align[queryIdx].emplace_back(alignInf, refIdx);
        }
    }

    for (auto &a : align) {
        std::sort(a.begin(), a.end());
    }
}

void Aligner::mergeAlignInf() {
    mergeAlignInfHelper(_readToContig, *_pReadToCtgDB, *_pReadDB, *_pContigDB);
    mergeAlignInfHelper(_readToRef, *_pReadToRefDB, *_pReadDB, *_pReferenceDB);
    mergeAlignInfHelper(_contigToRef, *_pCtgToRefDB, *_pContigDB, *_pReferenceDB);
}

void
Aligner::simpleAlign() {
    auto &contigDB = *_pContigDB;

    _alignReference.reset(contigDB);

    std::map<std::string, std::pair<std::size_t, std::size_t>> ctgCov;

    for (decltype(contigDB.size()) i = 0; i < contigDB.size(); ++i) {
        auto contigIdx = i;

        if (!_ctgFilterFlag[contigIdx]) { continue; }

        /*std::map<std::pair<std::size_t, bool>, std::vector<bool>> hitRecord;
        std::vector<std::vector<bool>> cov(2, std::vector<bool>(contigDB[contigIdx].size(), false));

        for (auto &align : _contigToRef[contigIdx]) {
            auto refIdx = align.getRefIdx();

            if (!_refFilterFlag[refIdx]) { continue; }

            auto ctgBegin = align.getBaseAlignInf().getQueryBegin();
            auto ctgEnd = align.getBaseAlignInf().getQueryEnd();
            auto forward = align.getBaseAlignInf().isForward();

            auto it = hitRecord.find({refIdx, forward});
            if (it == hitRecord.end()) {
                it = hitRecord.emplace(std::make_pair(refIdx, forward),
                                       std::vector<bool>(contigDB[contigIdx].size(), false)).first;
            }

            auto &arr = it->second;

            std::fill(arr.begin() + ctgBegin, arr.begin() + ctgEnd, true);
            //for (decltype(ctgBegin) j = ctgBegin; j < ctgEnd; ++j) {
            //    arr[j] = true;
            //}
        }

        std::set<std::pair<std::size_t, bool>> goodRefIdx;

        for (auto &kv : hitRecord) {
            auto refIdx = kv.first;
            auto &record = kv.second;

            auto total = record.size();
            auto cnt = std::accumulate(record.begin(), record.end(), static_cast<decltype(total)>(0));

            //std::cout << i << " " << cnt << std::endl;

            if (cnt * 1.0 / total >= _ctgToRefTotalRatio) {
                goodRefIdx.insert(refIdx);
            }
        }*/

        for (auto &align : _contigToRef[contigIdx]) {
            auto refIdx = align.getRefIdx();

            if (!_refFilterFlag[refIdx]) { continue; }

            auto forward = align.getBaseAlignInf().isForward();

            if (_ctgFilterForward[contigIdx] != forward) {
                continue;
            }

            //if (goodRefIdx.count({refIdx, forward}) == 0) { continue; }

            auto ctgBegin = align.getBaseAlignInf().getQueryBegin();
            auto ctgEnd = align.getBaseAlignInf().getQueryEnd();
            auto refBegin = align.getBaseAlignInf().getRefBegin();
            auto ctgLen = contigDB[contigIdx].size();
            auto queryDiff = align.getBaseAlignInf().getQueryDiff();
            auto refDiff = align.getBaseAlignInf().getRefDiff();

            //if ((ctgEnd - ctgBegin) * 1.0 / ctgLen < _ctgToRefRatio) { continue; }

            //if (ctgEnd - ctgBegin < _ctgToRefMinLen) { continue; }

            if (!forward) {
                flipPosition(ctgBegin, ctgEnd, ctgLen);
            }

            std::vector<pos_t> refPos;

            ParseAlignTools::exactAlign(ctgBegin, refBegin, true, queryDiff, refDiff, [&](std::size_t ctgCur, std::size_t refCur) {
                refPos.push_back(refCur);
            });

            _alignReference.insert(contigIdx, ctgBegin, ctgEnd, refIdx + 1, refPos, forward);

            // Calculate coverage
            //auto &arr = forward ? cov[0] : cov[1];
            //std::fill(arr.begin() + ctgBegin, arr.begin() + ctgEnd, true);
        }

        //ctgCov.emplace(std::piecewise_construct,
        //               std::forward_as_tuple(contigDB[contigIdx].name()),
        //               std::forward_as_tuple(std::accumulate(cov[0].begin(), cov[0].end(), static_cast<std::size_t>(0)),
        //                                     std::accumulate(cov[1].begin(), cov[1].end(), static_cast<std::size_t>(0)))
        //);

    }

    //return ctgCov;
}

void Aligner::addExtraPosition() {
    for (std::size_t i = 0; i < _ctgFilterFlag.size(); ++i) {
        if (_ctgFilterFlag[i]) {
            _alignReference.addExtraPosition(i, _ctgFilterForward[i]);
        }
    }
}

bool Aligner::setUseCtgAlign(const std::string &ctgName, bool forward, bool use) {
    auto id = _pContigDB->seqId(ctgName);
    if (id != _pContigDB->NOT_FOUND_ID) {
        _alignReference.setAlign(id, forward, use);
        return true;
    }

    return false;
}

bool Aligner::queryContig(std::size_t contigIdx, std::size_t contigPos, bool forward,
                          std::vector<std::pair<Aligner::ref_pos_t, Aligner::ref_pos_t>> &positions) const {
    auto query = _alignReference.query(contigIdx, (AlignReference::pos_t) contigPos, forward);

    for (auto &refPos : query) {
        positions.emplace_back(
                ref_pos_t(forward ? contigIdx + 1 : -contigIdx - 1, contigPos),
                refPos);
    }

    return !query.empty();
}

void Aligner::flipPosition(std::size_t &left, std::size_t &right, std::size_t length) {
    auto tmp = left;
    left = length - right;
    right = length - tmp;
}

Aligner::Aligner(std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
                 std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
                 std::shared_ptr<const ISeqDatabase<SeqInf>> pReferenceDB,
                 std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToContigDB,
                 std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToRefDB,
                 std::shared_ptr<const IAlignDatabase<AlignInf>> pContigToRefDB) :
        _pReadDB(std::move(pReadDB)),
        _pContigDB(std::move(pContigDB)),
        _pReferenceDB(std::move(pReferenceDB)),
        _pReadToCtgDB(std::move(pReadToContigDB)),
        _pReadToRefDB(std::move(pReadToRefDB)),
        _pCtgToRefDB(std::move(pContigToRefDB)),
        _isInit(false) {

    clearRefFilter(true);
    clearCtgFilter(true);
    _alignReference.reset(*_pContigDB);
    mergeAlignInf();
}

void Aligner::setReadToCtgTopK(int readToCtgTopK) {
    Aligner::_readToCtgTopK = readToCtgTopK;
}

void Aligner::setReadToRefTopK(int readToRefTopK) {
    Aligner::_readToRefTopK = readToRefTopK;
}

void Aligner::setCtgToRefTopK(int ctgToRefTopK) {
    Aligner::_ctgToRefTopK = ctgToRefTopK;
}

void Aligner::setReadToCtgRatio(double readToCtgRatio) {
    Aligner::_readToCtgRatio = readToCtgRatio;
}

void Aligner::setReadToRefRatio(double readToRefRatio) {
    Aligner::_readToRefRatio = readToRefRatio;
}

void Aligner::setCtgToRefRatio(double ctgToRefRatio) {
    Aligner::_ctgToRefRatio = ctgToRefRatio;
}

void Aligner::setCtgToRefTotalRatio(double ctgToRefTotalRatio) {
    Aligner::_ctgToRefTotalRatio = ctgToRefTotalRatio;
}

void Aligner::setCtgToRefMinLen(std::size_t ctgToRefMinLen) {
    Aligner::_ctgToRefMinLen = ctgToRefMinLen;
}

bool Aligner::setRefFilter(const std::string &refName, bool accepted) {
    auto id = _pReferenceDB->seqId(refName);
    if (id != _pReferenceDB->NOT_FOUND_ID) {
        _refFilterFlag[id] = accepted;
        return true;
    }

    return false;
}

void Aligner::clearRefFilter(bool accepted) {
    _refFilterFlag.assign(_pReferenceDB->size(), accepted);
}

bool Aligner::setCtgFilter(const std::string &refName, bool forward, bool accepted) {
    auto id = _pContigDB->seqId(refName);
    if (id != _pContigDB->NOT_FOUND_ID) {
        _ctgFilterFlag[id] = accepted;
        _ctgFilterForward[id] = forward;
        return true;
    }

    return false;
}

void Aligner::clearCtgFilter(bool accepted) {
    _ctgFilterFlag.assign(_pContigDB->size(), accepted);
    _ctgFilterForward.assign(_pContigDB->size(), true);
}
