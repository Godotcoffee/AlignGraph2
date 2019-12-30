//
// Created by amadeus on 19-5-6.
//

#ifndef PAGRAPH_ABSTRACTALIGNDATABASE_HPP
#define PAGRAPH_ABSTRACTALIGNDATABASE_HPP

#include <vector>
#include "IAlignDatabase.hpp"
#include "seq/ISeqDatabase.hpp"

class AbstractAlignDatabase : public IAlignDatabase<AlignInf> {

protected:
    std::vector<AlignInf> _alignments;

public:
    std::size_t size() const override;

    const AlignInf &operator[](std::size_t id) const override;

    AlignInf &operator[](std::size_t id) override;

    ~AbstractAlignDatabase() override = 0;
};


#endif //PAGRAPH_ABSTRACTALIGNDATABASE_HPP
