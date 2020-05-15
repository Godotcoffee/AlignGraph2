//
// Created by Godot on 2019/8/18.
//

#include <fstream>
#include <sstream>
#include <algorithm>
#include "align/AlignmentHelper.hpp"
#include "AlignData.hpp"

std::pair<std::size_t, std::size_t>
AlignData::sliceHelper(const std::string &tstr, std::size_t originStart, std::size_t sliceStart, std::size_t sliceEnd) {
    std::size_t left, right;
    std::size_t cnt = 0;

    for (left = 0; left < tstr.size(); ++left) {
        if (tstr[left] == '-') { continue; }
        if (originStart + cnt >= sliceStart) {
            break;
        }
        ++cnt;
    }

    for (right = left; right < tstr.size(); ++right) {
        if (tstr[right] == '-') { continue; }
        if (originStart + cnt >= sliceEnd) {
            break;
        }
        ++cnt;
    }

    return {left, right};
}

std::vector<std::vector<std::pair<dagcon::Alignment, unsigned>>>
AlignData::readFromRefFile(const std::string &filePath, const std::string &skeleton, std::size_t partLen) {
    if (skeleton.empty()) {
        return {};
    }
    auto partNum = (skeleton.size() + partLen - 1) / partLen;
    std::vector<std::vector<std::pair<dagcon::Alignment, unsigned>>> part(partNum);

    AlignmentHelper::loadFromRefFile(filePath, [&](const SimpleAlign &aln, const std::string &l1, const std::string &qstr, const std::string &tstr) {
        auto toStart = aln.refBegin;
        auto toEnd = aln.refEnd;
        auto &srcId = aln.queryName;
        auto &toId = aln.refName;
        auto forward = aln.forward;
        auto score = std::atoi(aln.score.c_str());

        auto leftPart = toStart / partLen;
        auto rightPart = std::min((toEnd - 1) / partLen, partNum - 1);

        for (auto i = leftPart; i <= rightPart; ++i) {
            dagcon::Alignment dagAln;
            dagAln.sid = srcId;
            dagAln.id = toId;
            dagAln.strand = forward ? '+' : '-';
            dagAln.start = i == leftPart ? toStart - leftPart * partLen + 1 : 1;
            dagAln.end = i == rightPart ? toEnd - rightPart * partLen + 1 : partLen;
            dagAln.tlen = partLen;

            auto slice = sliceHelper(tstr, toStart, i * partLen, std::min((i + 1) * partLen, skeleton.size()));

            dagAln.qstr = qstr.substr(slice.first, slice.second - slice.first);
            dagAln.tstr = tstr.substr(slice.first, slice.second - slice.first);

            auto newAln = normalizeGaps(dagAln);
            //trimAln(newAln);
            newAln.end = dagAln.end;

            part[i].emplace_back(newAln, score);
        }
    });

    return part;
}

std::vector<std::size_t> AlignData::weightAln(const std::vector<std::pair<dagcon::Alignment, unsigned>> &alignment, std::size_t alpha) {
    if (alignment.empty()) {
        return {};
    }

    std::vector<std::size_t> weights;

    auto minmax = std::minmax_element(alignment.begin(), alignment.end(),
                                      [](const std::pair<dagcon::Alignment, unsigned> &lhs,
                                         const std::pair<dagcon::Alignment, unsigned> &rhs) {
                                          return lhs.second < rhs.second;
                                      });

    auto maxScore = minmax.second->second;
    auto minScore = minmax.first ->second;
    auto scoreRange = maxScore - minScore;

    weights.reserve(alignment.size());
    for (auto &aln : alignment) {
        weights.push_back(
                std::max(static_cast<std::size_t>((aln.second - minScore) * 1.0 / std::max(scoreRange, 1U) * alpha), std::size_t(1))
        );
    }

    return weights;
}
