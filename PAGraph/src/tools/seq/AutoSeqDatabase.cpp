//
// Created by amadeus on 19-12-26.
//

#include "AutoSeqDatabase.hpp"
#include "seq/SeqHelper.hpp"

AutoSeqDatabase::AutoSeqDatabase(const std::string &filePath) {
    std::size_t count = 0;
    SeqHelper::autoLoadFromFile(filePath, [&](const std::string &comment, const std::string &seq) {
        std::string name = comment.substr(1);
        _seqs.emplace_back(seq, name);
        _nameToId[name] = count++;
    });
}
