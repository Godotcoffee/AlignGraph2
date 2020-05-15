//
// Created by Godot on 2019/8/18.
//

#ifndef PBDAGTEST_ALIGNDATA_HPP
#define PBDAGTEST_ALIGNDATA_HPP


#include <vector>
#include "Alignment.hpp"

class AlignData {
private:
    static std::pair<std::size_t, std::size_t>
    sliceHelper(const std::string &tstr, std::size_t originStart, std::size_t sliceStart, std::size_t sliceEnd);

public:
    static std::vector<std::vector<std::pair<dagcon::Alignment, unsigned>>>
    readFromRefFile(const std::string &filePath, const std::string &skeleton, std::size_t partLen);

    static std::vector<std::size_t> weightAln(const std::vector<std::pair<dagcon::Alignment, unsigned>> &alignment, std::size_t alpha = 250);
};


#endif //PBDAGTEST_ALIGNDATA_HPP
