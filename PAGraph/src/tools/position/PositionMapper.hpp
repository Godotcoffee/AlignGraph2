//
// Created by amadeus on 19-5-27.
//

#ifndef PAGRAPH_POSITIONMAPPER_HPP
#define PAGRAPH_POSITIONMAPPER_HPP


#include <seq/ISeqDatabase.hpp>
#include <seq/SeqInf.hpp>

class PositionMapper {
private:
    std::vector<std::size_t> _sizes;
    std::vector<std::size_t> _startPos;

    static void generateStartPosHelper(std::vector<std::size_t> &start, const ISeqDatabase<SeqInf> &seqDB);
public:
    explicit PositionMapper(const ISeqDatabase<SeqInf> &ref);

    std::size_t extraStart() const;

    std::size_t dualToSingle(std::int64_t refIdx, std::int64_t pos) const;

    std::pair<std::int64_t, std::int64_t> singleToDual(std::size_t singlePos) const;

    std::size_t size(std::int64_t refIdx) const;
};


#endif //PAGRAPH_POSITIONMAPPER_HPP
