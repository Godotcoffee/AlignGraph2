//
// Created by amadeus on 19-5-10.
//

#ifndef PAGRAPH_ABRUIJNNODE_HPP
#define PAGRAPH_ABRUIJNNODE_HPP

#include "KMerAdjNode.hpp"

template<typename PosType, typename IndexType, typename CountType = unsigned short>
class ABruijnNode {
private:
    typedef KMerAdjNode<PosType, IndexType, CountType> Node;
    Node *_node;            // Use pointer instead of reference for copy.
    std::size_t _indexVal;
public:
    explicit ABruijnNode() : _node(nullptr), _indexVal(static_cast<std::size_t>(-1)) {}

    explicit ABruijnNode(Node &node, const std::size_t &index) : _node(&node), _indexVal(index) {}

    std::size_t size() const {
        return _node->getAllPositions().size();
    }

    const PosType &operator[](std::size_t i) const {
        return _node->getAllPositions()[i];
    }

    const CountType getPosAbundance(std::size_t index) const {
        return _node->getPosCount(index);
    }

    void setUsed(std::size_t index) {
        _node->setFlag(index);
    }

    void unsetUsed(std::size_t index) {
        _node->unsetFlag(index);
    }

    bool isUsed(std::size_t index) const {
        return _node->getFlag(index);
    }

    void setIndex(const std::size_t index) {
        _indexVal = index;
    }

    std::size_t getIndex() const {
        return _indexVal;
    }
};


#endif //PAGRAPH_ABRUIJNNODE_HPP
