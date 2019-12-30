//
// Created by amadeus on 19-5-6.
//

#include "AbstractAlignDatabase.hpp"

std::size_t AbstractAlignDatabase::size() const {
    return _alignments.size();
}

const AlignInf &AbstractAlignDatabase::operator[](std::size_t id) const {
    return _alignments[id];
}

AlignInf &AbstractAlignDatabase::operator[](std::size_t id) {
    return
            const_cast<AlignInf &>(
                    static_cast<const AbstractAlignDatabase &>(*this)[id]
            );
}

AbstractAlignDatabase::~AbstractAlignDatabase() = default;
