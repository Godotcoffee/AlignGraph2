//
// Created by amadeus on 19-6-24.
//

#include "ParseAlignTools.hpp"

void
ParseAlignTools::parseDiff(const std::string &queryDiffStr, const std::string &refDiffStr,
                           std::vector<bool> &queryDiff, std::vector<bool> &refDiff) {

    for (std::size_t i = 0; i < queryDiffStr.size(); ++i) {
        if (queryDiffStr[i] == '-') {
            queryDiff.push_back(true);
            refDiff.push_back(false);
        } else if (refDiffStr[i] == '-') {
            queryDiff.push_back(false);
            refDiff.push_back(true);
        } else if (queryDiffStr[i] != refDiffStr[i]) {
            queryDiff.push_back(true);
            refDiff.push_back(true);
        } else {
            queryDiff.push_back(false);
            refDiff.push_back(false);
        }
    }
}
