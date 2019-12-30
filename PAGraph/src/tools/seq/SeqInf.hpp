//
// Created by amadeus on 19-5-6.
//

#ifndef PAGRAPH_SEQINF_HPP
#define PAGRAPH_SEQINF_HPP

#include <cstddef>
#include <string>
#include "CompressedSeq.hpp"

class SeqInf {
private:
    CompressedSeq _compressedSeq;
    std::string _comment;
public:
    explicit SeqInf();

    explicit SeqInf(const std::string &seqStr, std::string comment);

    std::size_t size() const;

    std::string toString(bool forward = true) const;

    const std::string &getComment() const;

    void setComment(const std::string &comment);

    char quickBaseAt(std::size_t idx, bool forward = true) const;

    const std::string &name() const;
};


#endif //PAGRAPH_SEQINF_HPP
