//
// Created by amadeus on 19-5-6.
//

#ifndef PAGRAPH_ABSTRACTSEQDATABASE_HPP
#define PAGRAPH_ABSTRACTSEQDATABASE_HPP

#include <cstddef>
#include <unordered_map>

#include "ISeqDatabase.hpp"
#include "SeqInf.hpp"

/**
 * Abstract class of ISeqDatabase. Subclass should add data to _seqs and _nameToId
 */
class AbstractSeqDatabase : public ISeqDatabase<SeqInf> {

protected:
    std::vector<SeqInf> _seqs;
    std::unordered_map<std::string, std::size_t> _nameToId;

public:
    std::size_t size() const override;

    const SeqInf &operator[](std::size_t id) const override;

    SeqInf &operator[](std::size_t id) override;

    std::string seqName(std::size_t id) const override;

    bool containsName(const std::string &name) const override;

    std::size_t seqId(const std::string &name) const override;

    ~AbstractSeqDatabase() override = 0;
};


#endif //PAGRAPH_ABSTRACTSEQDATABASE_HPP
