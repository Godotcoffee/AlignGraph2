//
// Created by amadeus on 19-5-10.
//

#ifndef PAGRAPH_PABRUIJNNODE_HPP
#define PAGRAPH_PABRUIJNNODE_HPP

#include <cstddef>
#include <cstdint>
#include "ABruijnNode.hpp"

template<typename PosType, typename IndexType, typename CountType = unsigned short>
class PABruijnNode {
public:
    typedef std::pair<std::size_t, std::size_t> UniqueType;
private:
    typedef ABruijnNode<PosType, IndexType> Node;
    Node _node;
    std::size_t _posIndex;
public:
    explicit PABruijnNode() : _posIndex(static_cast<std::size_t>(-1)) {}

    explicit PABruijnNode(const Node &node, std::size_t posIndex) : _node(node), _posIndex(posIndex) {}

    Node &getABruijnNode() {
        return _node;
    }

    const Node &getABruijnNode() const {
        return _node;
    }

    void setPosIndex(std::size_t pos) {
        _posIndex = pos;
    }

    std::size_t getPosIndex() const {
        return _posIndex;
    }

    bool isExist() const {
        return _posIndex < _node.size();
    }

    PosType getPosition() const {
        return isExist() ? _node[_posIndex] : PosType();
    }

    CountType getAbundance() const {
        return isExist() ? _node.getPosAbundance(_posIndex) : CountType();
    }

    void setUsed() {
        if (isExist()) {
            _node.setUsed(_posIndex);
        }
    }

    void unsetUsed() {
        if (isExist()) {
            _node.unsetUsed(_posIndex);
        }
    }

    bool isUsed() const {
        return isExist() ? _node.isUsed(_posIndex) : true;
    }

    UniqueType getUniqueHelper() const {
        return {getABruijnNode().getIndex(), _posIndex};
    }
};


#endif //PAGRAPH_PABRUIJNNODE_HPP
