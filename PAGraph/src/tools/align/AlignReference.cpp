//
// Created by amadeus on 19-4-18.
//

#include "AlignReference.hpp"

void AlignReference::reset(const ISeqDatabase<SeqInf> &queryDB) {
    _fPositions.clear();
    _rPositions.clear();
    for (std::size_t i = 0; i < queryDB.size(); ++i) {
        _fPositions.emplace_back(queryDB[i].size());
        _rPositions.emplace_back(queryDB[i].size());
    }

    _fIsAlign.assign(queryDB.size(), false);
    _rIsAlign.assign(queryDB.size(), false);
}

bool
AlignReference::insert(std::size_t queryIdx, AlignReference::pos_t start, AlignReference::pos_t end, std::size_t refIdx,
                       AlignReference::pos_t posStart, AlignReference::pos_t posEnd, bool forward) {
    auto &positions = forward ? _fPositions[queryIdx] : _rPositions[queryIdx];
    auto left = start;
    auto right = end;

    //double cur = posStart;
    //double step = (posEnd - posStart) * 1.0 / (end - start);
    auto cur = posStart;

    for (auto i = left; i < right; ++i) {
        //positions[i].emplace_back(refIdx, (pos_t) std::round(cur));
        //cur += step;
        positions[i].emplace_back(refIdx, cur);
        ++cur;
    }
    (forward ? _fIsAlign[queryIdx] : _rIsAlign[queryIdx]) = true;
    return true;
}

bool
AlignReference::insert(std::size_t queryIdx, AlignReference::pos_t start, AlignReference::pos_t end, std::size_t refIdx,
                       const std::vector<AlignReference::pos_t> &refPos, bool forward) {
    auto &positions = forward ? _fPositions[queryIdx] : _rPositions[queryIdx];
    auto left = start;
    auto right = end;

    //double cur = posStart;
    //double step = (posEnd - posStart) * 1.0 / (end - start);

    for (auto i = left; i < right; ++i) {
        //positions[i].emplace_back(refIdx, (pos_t) std::round(cur));
        //cur += step;
        positions[i].emplace_back(refIdx, refPos[i - left]);
    }
    (forward ? _fIsAlign[queryIdx] : _rIsAlign[queryIdx]) = true;
    return true;
}

const std::vector<AlignReference::ref_pos_t> &
AlignReference::query(std::size_t queryIdx, AlignReference::pos_t pos, bool forward) const {
    auto &positions = forward ? _fPositions[queryIdx] : _rPositions[queryIdx];
    if (pos >= 0 && static_cast<std::size_t>(pos) < positions.size()) {
        return positions[pos];
    }

    return _empty_pos;
}

void AlignReference::addExtraPosition(std::size_t queryIdx, bool forward) {
    auto &pos = forward ? _fPositions[queryIdx] : _rPositions[queryIdx];

    for (auto &p : pos) {
        if (p.empty()) {
            p.emplace_back(std::int64_t(0), pos_t(0));
        }
    }


}

std::size_t AlignReference::size() const {
    return _fPositions.size();
}

std::size_t AlignReference::size(std::size_t idx) const {
    return _fPositions[idx].size();
}

bool AlignReference::isAlign(std::size_t refIdx, bool forward) {
    return forward ? _fIsAlign[refIdx] : _rIsAlign[refIdx];
}

void AlignReference::setAlign(std::size_t refIdx, bool forward, bool align) {
    (forward ? _fIsAlign[refIdx] : _rIsAlign[refIdx]) = align;
}