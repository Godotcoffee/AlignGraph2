//
// Created by amadeus on 19-5-9.
//

#ifndef PAGRAPH_ABSTRACTKMERDATABASE_HPP
#define PAGRAPH_ABSTRACTKMERDATABASE_HPP

#include <vector>
#include "IKmerDatabase.hpp"
#include "KmerInf.hpp"

class AbstractKmerDatabase : public IKmerDatabase<KmerInf> {

protected:
    std::vector<KmerInf> _kmers;

    std::size_t _kmerSize = 0;

public:
    size_t size() const override;

    const KmerInf &operator[](std::size_t id) const override;

    KmerInf &operator[](std::size_t id) override;

    size_t kmerSize() const override;

    ~AbstractKmerDatabase() override = 0;
};


#endif //PAGRAPH_ABSTRACTKMERDATABASE_HPP
