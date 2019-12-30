//
// Created by amadeus on 19-5-6.
//

#include "AbstractSeqDatabase.hpp"
#include <limits>

std::size_t AbstractSeqDatabase::size() const {
    return _seqs.size();
}

const SeqInf &AbstractSeqDatabase::operator[](std::size_t id) const {
    return _seqs[id];
}

SeqInf &AbstractSeqDatabase::operator[](std::size_t id) {
    return _seqs[id];
}

std::string AbstractSeqDatabase::seqName(std::size_t id) const {
    return _seqs[id].name();
}

bool AbstractSeqDatabase::containsName(const std::string &name) const {
    return _nameToId.count(name) > 0;
}

std::size_t AbstractSeqDatabase::seqId(const std::string &name) const {
    auto it = _nameToId.find(name);
    return it != _nameToId.end() ? it->second : NOT_FOUND_ID;
}

AbstractSeqDatabase::~AbstractSeqDatabase() = default;
