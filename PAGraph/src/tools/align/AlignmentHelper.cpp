//
// Created by amadeus on 19-11-12.
//

#include "AlignmentHelper.hpp"
#include <sstream>
#include <fstream>
#include <vector>

void
AlignmentHelper::loadFromRefFile(const std::string &refPath,
                     std::function<void(SimpleAlign &, const std::string &line1, const std::string &line2,
                                        const std::string &line3)> f) {
    std::ifstream inAlign(refPath);
    if (inAlign.is_open()) {
        std::string line;
        std::stringstream ss;

        std::string queryName, refName, forward, score;
        std::size_t queryBegin = 0, queryEnd = 0, querySize = 0, refBegin = 0, refEnd = 0, refSize = 0;
        SimpleAlign simpleAlign {"", "", false, "", 0, 0, 0, 0, 0, 0};
        std::string line1, line2, line3;

        for (std::size_t lineCount = 0; ; ++lineCount) {
            if (lineCount % 3 == 0) {
                if (!std::getline(inAlign, line1)) { break; }
                ss.clear();
                ss.str(line1);
                ss >> queryName >> refName >> forward >> score >> queryBegin >> queryEnd >> querySize >> refBegin
                   >> refEnd >> refSize;
                if (!ss.fail()) {
                    simpleAlign = {
                            queryName, refName, forward == "F", score, queryBegin, queryEnd, querySize, refBegin, refEnd, refSize
                    };
                } else {
                    simpleAlign = {"", "", false, "", 0, 0, 0, 0, 0, 0};
                }
            } else if (lineCount % 3 == 1) {
                if (!std::getline(inAlign, line2)) { break; }
            } else {
                if (!std::getline(inAlign, line3)) { break; }
                f(simpleAlign, line1, line2, line3);
            }
        }
        inAlign.close();
    }
}

void
AlignmentHelper::simpleLoadFromRefFile(const std::string &refPath, std::function<void(SimpleAlign &)> f) {
    loadFromRefFile(refPath, [&](SimpleAlign &align, const std::string &, const std::string &, const std::string &) {
        f(align);
    });
}

void
AlignmentHelper::filterRefFile(const std::string &refPath, const std::string &outPath, std::function<bool(SimpleAlign&)> f) {

    std::ofstream outAlign(outPath);
    if (outAlign.is_open()) {
        loadFromRefFile(refPath, [&](SimpleAlign &align, const std::string &l1, const std::string &l2, const std::string &l3) {
            if (f(align)) {
                outAlign << l1 << std::endl;
                outAlign << l2 << std::endl;
                outAlign << l3 << std::endl;
            }
        });
        outAlign.close();
    }
}
