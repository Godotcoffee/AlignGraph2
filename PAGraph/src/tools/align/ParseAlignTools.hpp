//
// Created by amadeus on 19-6-24.
//

#ifndef PAGRAPH_PARSEALIGNTOOLS_HPP
#define PAGRAPH_PARSEALIGNTOOLS_HPP

#include <vector>
#include <string>

class ParseAlignTools {
public:
    static void parseDiff(const std::string &queryDiffStr, const std::string &refDiffStr,
                          std::vector<bool> &queryDiff, std::vector<bool> &refDiff);

    template <typename T>
    static void exactAlign(
            std::size_t queryBegin, std::size_t refBegin, bool forward,
            const std::vector<bool> &queryDiff, const std::vector<bool> &refDiff, T functor);
};

#include "ParseAlignTools.tcc"


#endif //PAGRAPH_PARSEALIGNTOOLS_HPP
