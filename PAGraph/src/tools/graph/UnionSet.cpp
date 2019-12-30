//
// Created by amadeus on 19-12-5.
//

#include "UnionSet.hpp"

void UnionSet::init(std::size_t size) {
    _parent.resize(size);
    for (std::size_t i = 0; i < size; ++i) {
        _parent[i] = i;
    }
}

std::size_t UnionSet::find(std::size_t pos) {
    return pos == _parent[pos] ? pos : (_parent[pos] = find(_parent[pos]));
}

void UnionSet::unionTwo(std::size_t left, std::size_t right) {
    _parent[find(right)] = find(left);
}
