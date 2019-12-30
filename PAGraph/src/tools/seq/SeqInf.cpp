//
// Created by amadeus on 19-7-5.
//

#include "SeqInf.hpp"

SeqInf::SeqInf()
        : _compressedSeq(""), _comment("") {}

SeqInf::SeqInf(const std::string &seqStr, std::string comment)
        : _compressedSeq(seqStr), _comment(std::move(comment)) {}

std::size_t SeqInf::size() const {
    return _compressedSeq.size();
}

std::string SeqInf::toString(bool forward) const {
    return _compressedSeq.toString(forward);
}

const std::string &SeqInf::getComment() const {
    return _comment;
}

void SeqInf::setComment(const std::string &comment) {
    SeqInf::_comment = comment;
}

char SeqInf::quickBaseAt(std::size_t idx, bool forward) const {
    return _compressedSeq.baseAt(idx, forward);
}

const std::string &SeqInf::name() const {
    return getComment();
}

