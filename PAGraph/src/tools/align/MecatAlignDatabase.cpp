//
// Created by amadeus on 19-5-5.
//

#include "MecatAlignDatabase.hpp"
#include "AlignmentHelper.hpp"

MecatAlignDatabase::MecatAlignDatabase(const std::string &filePath) {
    AlignmentHelper::loadFromRefFile(filePath, [&](SimpleAlign &align, const std::string &l1,
            const std::string &l2, const std::string &l3) {
        std::vector<bool> readDiff, refDiff;

        ParseAlignTools::parseDiff(l2, l3, readDiff, refDiff);

        _alignments.emplace_back(align.queryName, align.refName, std::atoll(align.score.c_str()), align.queryBegin, align.queryEnd,
                align.refBegin, align.refEnd, align.forward, readDiff, refDiff);
    });

    std::sort(_alignments.begin(), _alignments.end());
}