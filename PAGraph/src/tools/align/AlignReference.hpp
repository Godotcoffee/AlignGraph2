//
// Created by amadeus on 19-4-18.
//

#ifndef PAGRAPH_ALIGNREFERENCE_HPP
#define PAGRAPH_ALIGNREFERENCE_HPP

#include <vector>
#include <cmath>
#include <cstdint>
#include <seq/SeqInf.hpp>
#include <seq/ISeqDatabase.hpp>
#include <iostream>

class AlignReference {
public:
    typedef std::int64_t pos_t;
    typedef std::pair<std::int64_t, pos_t> ref_pos_t;
private:
    std::vector<std::vector<std::vector<ref_pos_t>>> _fPositions;
    std::vector<std::vector<std::vector<ref_pos_t>>> _rPositions;

    std::vector<bool> _fIsAlign;
    std::vector<bool> _rIsAlign;

    std::vector<AlignReference::ref_pos_t> _empty_pos;
public:
    void reset(const ISeqDatabase<SeqInf> &queryDB);

    bool insert(std::size_t queryIdx, pos_t start, pos_t end, std::size_t refIdx, pos_t posStart, pos_t posEnd,
                bool forward = true);

    bool insert(std::size_t queryIdx, pos_t start, pos_t end, std::size_t refIdx, const std::vector<pos_t> &refPos,
                bool forward = true);

    const std::vector<ref_pos_t> &query(std::size_t queryIdx, pos_t pos, bool forward = true) const;

    void addExtraPosition(std::size_t queryIdx, bool forward);

    std::size_t size() const;

    std::size_t size(std::size_t idx) const;

    bool isAlign(std::size_t refIdx, bool forward = true);

    void setAlign(std::size_t refIdx, bool forward, bool align);
};


#endif //GATB_MYSELF_ALIGNREFERENCE_HPP
