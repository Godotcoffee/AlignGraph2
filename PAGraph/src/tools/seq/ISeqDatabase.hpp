//
// Created by amadeus on 19-5-5.
//

#ifndef PAGRAPH_ISEQDATABASE_HPP
#define PAGRAPH_ISEQDATABASE_HPP

#include <cstddef>
#include <string>
#include <limits>
#include "thread/IMultiThread.hpp"

/**
 * Interface of a set of Sequences.
 */
template <typename T>
class ISeqDatabase : public IMultiThread<T> {
public:
    /**
     * Use special number to indicate NOT FOUND
     */
    static const auto NOT_FOUND_ID = std::numeric_limits<std::size_t>::max();

    /**
     * @return Number of sequences.
     */
    virtual std::size_t size() const = 0;

    /**
     * Const getter
     * @return Sequence information with a specific id.
     */
    virtual const T& operator[](std::size_t id) const = 0;

    /**
     * Getter
     * @return Sequence information with a specific id.
     */
    virtual T& operator[](std::size_t id) = 0;

    /**
     * @return Name of a specific sequence. Empty if doesn't exist.
     */
    virtual std::string seqName(std::size_t id) const = 0;

    /**
     * @return Whether contains the sequence with a specific name.
     */
    virtual bool containsName(const std::string &name) const = 0;

    /**
     * @return Sequence id with a specific name. Max value of size_t if doesn't exist.
     */
    virtual std::size_t seqId(const std::string &name) const = 0;

    virtual ~ISeqDatabase() = default;
};


#endif //PAGRAPH_ISEQDATABASE_HPP
