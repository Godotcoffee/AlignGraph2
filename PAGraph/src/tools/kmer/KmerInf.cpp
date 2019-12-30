//
// Created by amadeus on 19-7-5.
//

#include "KmerInf.hpp"

KmerInf::KmerInf(const std::string &kmer, decltype(KmerInf::_abundance) abundance) : _kmerSeq(kmer), _abundance(abundance) {}

std::string KmerInf::toString() const {
    return _kmerSeq.toString();
}

void KmerInf::setAbundance(decltype(KmerInf::_abundance) abundance) {
    _abundance = abundance;
}

decltype(KmerInf::_abundance) KmerInf::getAbundance() const {
    return _abundance;
}
