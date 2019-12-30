//
// Created by amadeus on 19-5-27.
//

#include "PositionMapper.hpp"
#include <algorithm>

PositionMapper::PositionMapper(const ISeqDatabase<SeqInf> &ref) {
    generateStartPosHelper(_startPos, ref);

    for (decltype(ref.size()) i = 0; i < ref.size(); ++i) {
        _sizes.push_back(ref[i].size());
    }
}

void PositionMapper::generateStartPosHelper(std::vector<std::size_t> &start, const ISeqDatabase<SeqInf> &seqDB) {
    start.clear();

    if (seqDB.size() == 0) {
        return;
    }

    start.push_back(seqDB[0].size());

    for (decltype(seqDB.size()) i = 1; i < seqDB.size(); ++i) {
        start.push_back(start.back() + 3 * seqDB[i - 1].size() +
                        std::max(seqDB[i - 1].size(), seqDB[i].size()));
    }

    start.push_back(start.back() + 4 * seqDB[seqDB.size() - 1].size());
}

std::size_t PositionMapper::extraStart() const {
    return !_startPos.empty() ? _startPos.back() : 0;
}

std::size_t PositionMapper::dualToSingle(std::int64_t refIdx, std::int64_t pos) const {
    if (refIdx == 0) { return 0; }
    auto idx = refIdx > 0 ? refIdx - 1 : -refIdx - 1;
    auto offset = refIdx > 0 ? 0 : 2 * _sizes[idx];
    return _startPos[idx] + offset + pos;
}

std::pair<std::int64_t, std::int64_t> PositionMapper::singleToDual(std::size_t singlePos) const {
    if (singlePos == 0) {
        return {0, 0};
    }
    auto it = std::upper_bound(_startPos.begin(), _startPos.end(), singlePos);
    if (it != _startPos.begin()) {
        it = std::prev(it);
    }

    auto idx = it - _startPos.begin();
    auto offset = singlePos - *it;

    if (offset >= static_cast<std::int64_t>(2 * _sizes[idx])) {
        offset -= static_cast<std::int64_t>(2 * _sizes[idx]);
        idx = -(idx + 1);
    } else {
        ++idx;
    }

    return {idx, offset};
}

std::size_t PositionMapper::size(std::int64_t refIdx) const {
    if (refIdx == 0) { return 0; }
    auto idx = refIdx > 0 ? refIdx - 1 : -refIdx - 1;
    return _sizes[idx];
}
