//
// Created by amadeus on 19-12-26.
//

#include <sstream>
#include "AutoSeqDatabase.hpp"
#include "seq/SeqHelper.hpp"

AutoSeqDatabase::AutoSeqDatabase(const std::string &filePath) {
    std::size_t count = 0;
    std::stringstream ss;
    std::string name;
    SeqHelper::autoLoadFromFile(filePath, [&](const std::string &comment, const std::string &seq) {
        ss.str("");
        ss.clear();
        ss << comment;
        ss >> name;
        std::string sp = name.substr(1);
        _seqs.emplace_back(seq, sp);
        _nameToId[sp] = count++;
    });
}
