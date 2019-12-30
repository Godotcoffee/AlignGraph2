//
// Created by amadeus on 19-6-24.
//

#include "MummerAlignDatabaseV2.hpp"

MummerAlignDatabaseV2::MummerAlignDatabaseV2(const std::string &filePath) {
    std::ifstream inalign(filePath);
    if (inalign.is_open()) {
        std::string line;
        std::stringstream ss;

        std::string queryName, refName, forward, readDiffStr, ignored;
        std::size_t queryBegin = 0, queryEnd = 0, refBegin = 0, refEnd = 0;
        bool bad = false;

        for (std::size_t lineCount = 0; std::getline(inalign, line); ++lineCount) {
            if (lineCount % 3 == 0) {
                ss.clear();
                ss.str(line);
                ss >> queryName >> refName >> forward >> ignored >> queryBegin >> queryEnd >> ignored >> refBegin
                   >> refEnd;
                if (ss.fail()) {
                    bad = true;
                }
            } else if (lineCount % 3 == 1) {
                if (!bad) {
                    readDiffStr = line;
                }
            } else {
                if (!bad) {
                    std::string &refDiffStr = line;

                    std::vector<bool> readDiff, refDiff;

                    ParseAlignTools::parseDiff(readDiffStr, refDiffStr, readDiff, refDiff);

                    _alignments.emplace_back(queryName, refName, queryEnd - queryBegin, queryBegin, queryEnd,
                                             refBegin, refEnd, forward == "F", readDiff, refDiff);
                }

                bad = false;
            }
        }
        inalign.close();
    }

    std::sort(_alignments.begin(), _alignments.end());
}
