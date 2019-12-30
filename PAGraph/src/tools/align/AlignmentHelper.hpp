//
// Created by amadeus on 19-11-12.
//

#ifndef PAGRAPH_ALIGNMENTHELPER_HPP
#define PAGRAPH_ALIGNMENTHELPER_HPP

#include <string>
#include <functional>

struct SimpleAlign {
    std::string queryName;
    std::string refName;
    bool forward;
    std::string score;
    std::size_t queryBegin;
    std::size_t queryEnd;
    std::size_t querySize;
    std::size_t refBegin;
    std::size_t refEnd;
    std::size_t refSize;
};

class AlignmentHelper {
public:
    static void simpleLoadFromRefFile(const std::string &refPath, std::function<void(SimpleAlign &)> f);

    static void filterRefFile(const std::string &refPath, const std::string &outPath, std::function<bool(SimpleAlign&)> f);

    static void loadFromRefFile(const std::string &refPath,
                                std::function<void(SimpleAlign &, const std::string &line1, const std::string &line2,
                                                   const std::string &line3)> f);
};


#endif //PAGRAPH_ALIGNMENTHELPER_HPP
