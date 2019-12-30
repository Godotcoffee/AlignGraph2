//
// Created by amadeus on 19-5-9.
//

#ifndef PAGRAPH_KMERINF_HPP
#define PAGRAPH_KMERINF_HPP

#include <seq/CompressedSeq.hpp>

class KmerInf {
private:
    CompressedSeq _kmerSeq;
    unsigned _abundance;
public:
    explicit KmerInf(const std::string &kmer, decltype(_abundance) abundance);

    std::string toString() const;

    void setAbundance(decltype(_abundance) abundance);

    decltype(_abundance) getAbundance() const;
};


#endif //PAGRAPH_KMERINF_HPP
