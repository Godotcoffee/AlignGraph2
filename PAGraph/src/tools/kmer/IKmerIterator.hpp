//
// Created by amadeus on 19-5-13.
//

#ifndef PAGRAPH_IKMERITERATOR_HPP
#define PAGRAPH_IKMERITERATOR_HPP

#include <functional>
#include <string>

class IKmerIterator {
public:
    virtual std::size_t kSize() const = 0;

    virtual void iterate(std::function<void(const std::uint64_t)> f, unsigned threadNum) const = 0;

    virtual ~IKmerIterator() = default;
};


#endif //PAGRAPH_IKMERITERATOR_HPP
