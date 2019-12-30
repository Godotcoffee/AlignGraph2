//
// Created by amadeus on 19-5-9.
//

#include "AbstractKmerDatabase.hpp"

size_t AbstractKmerDatabase::size() const {
    return _kmers.size();
}

const KmerInf &AbstractKmerDatabase::operator[](std::size_t id) const {
    return _kmers[id];
}

KmerInf &AbstractKmerDatabase::operator[](std::size_t id) {
    return
            const_cast<KmerInf &>(
                    static_cast<const AbstractKmerDatabase &>(*this)[id]
            );
}

size_t AbstractKmerDatabase::kmerSize() const {
    return _kmerSize;
}

AbstractKmerDatabase::~AbstractKmerDatabase() = default;
