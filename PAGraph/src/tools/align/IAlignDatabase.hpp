//
// Created by amadeus on 19-5-5.
//

#ifndef PAGRAPH_IALIGNDATABASE_HPP
#define PAGRAPH_IALIGNDATABASE_HPP

#include "AlignInf.hpp"

template <typename T>
class IAlignDatabase {
public:
    /**
     * @return Number of alignments
     */
    virtual std::size_t size() const = 0;

    /**
     * const getter
     * @return Alignment information with a specific id.
     */
    virtual const T &operator[](std::size_t id) const = 0;

    /**
     * getter
     * @return Alignment information with a specific id.
     */
    virtual T &operator[](std::size_t id) = 0;

    virtual ~IAlignDatabase() = default;
};


#endif //PAGRAPH_IALIGNDATABASE_HPP
