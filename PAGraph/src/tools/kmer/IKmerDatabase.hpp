//
// Created by amadeus on 19-5-9.
//

#ifndef PAGRAPH_IKMERDATABASE_HPP
#define PAGRAPH_IKMERDATABASE_HPP

#include <cstddef>
#include "thread/IMultiThread.hpp"

template<typename T>
class IKmerDatabase : public IMultiThread<T> {
public:
    /**
     * @return Number of k-mers
     */
    virtual std::size_t size() const = 0;

    /**
     * @return Kmer information
     */
    virtual const T &operator[](std::size_t id) const = 0;

    /**
     * @return Kmer information
     */
    virtual T &operator[](std::size_t id) = 0;

    /**
     * @return Value of K
     */
    virtual std::size_t kmerSize() const = 0;

    virtual ~IKmerDatabase() = default;
};


#endif //PAGRAPH_IKMERDATABASE_HPP
