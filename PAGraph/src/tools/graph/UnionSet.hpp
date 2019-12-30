//
// Created by amadeus on 19-12-5.
//

#ifndef PAGRAPH_UNIONSET_HPP
#define PAGRAPH_UNIONSET_HPP

#include <vector>

class UnionSet {
private:
    std::vector<std::size_t> _parent;
public:
    void init(std::size_t size);

    std::size_t find(std::size_t pos);

    void unionTwo(std::size_t left, std::size_t right);
};


#endif //PAGRAPH_UNIONSET_HPP
