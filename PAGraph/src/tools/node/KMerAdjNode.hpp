//
// Created by amadeus on 19-5-9.
//

#ifndef PAGRAPH_KMERADJNODE_HPP
#define PAGRAPH_KMERADJNODE_HPP

#include <iostream>
#include <vector>
#include <mutex>
#include <algorithm>
#include <iterator>
#include <functional>

template<typename PosType, typename IndexType, typename CountType = unsigned short>
class KMerAdjNode {
private:
    //const static unsigned _bitLen = sizeof(BitType) * 8;
    std::vector<IndexType> _children;
    std::vector<PosType> _positions;
    std::vector<CountType> _posCount;
    std::vector<bool> _usedFlag;
    mutable std::mutex _mutex;

    template<typename T>
    static void addWithMutex(std::mutex &mutex, std::vector<T> &store, const T &item);

    template<typename T>
    static void addDirectly(std::vector<T> &store, const T &item);

    template<typename T>
    static void addItem(std::mutex &mutex, std::vector<T> &store, const T &item, bool useMutex = false);

    template<typename T>
    static void addItems(std::mutex &mutex, std::vector<T> &store, const std::vector<T> &item, bool useMutex = false);

    template <typename T>
    static std::size_t removeDuplicate(std::vector<T> &array, std::function<bool(const T&, const T&)> cmp, bool shrink = false);

    template <typename T1, typename T2>
    static std::size_t cluster(
            std::vector<T1> &array, std::vector<T2> &counts, std::function<bool(const T1&, const T1&)> cmp, bool shrink = false);

    template <typename T1, typename T2>
    static void sortWithCount(std::vector<T1> &array, std::vector<T2> &count);

    void reset();
public:
    explicit KMerAdjNode();

    void setAbundance(int abundance);

    void addChild(const IndexType &child, bool useMutex = true);

    std::size_t mergeChild(std::function<bool(const IndexType&, const IndexType&)> cmp, bool shrink = false);

    void addPosition(const PosType &pos, bool useMutex = true);

    void addPosition(const std::vector<PosType> &pos, bool useMutex = true);

    std::size_t mergePosition(std::function<bool(const PosType&, const PosType&)> cmp, bool shrink = false);

    void sortPosition();

    CountType getPosCount(std::size_t index, bool useMutex = false);

    void resetUsedFlag(bool useMutex = false);

    void setFlag(std::size_t index, bool useMutex = false);

    void unsetFlag(std::size_t index, bool useMutex = false);

    bool getFlag(std::size_t index, bool useMutex = false) const;

    void resetAll(bool useMutex = false);

    const std::vector<IndexType> &getAllChild() const;

    const std::vector<PosType> &getAllPositions() const;

    const std::vector<CountType> &getAllCount() const;


};

#include "KMerAdjNode.tcc"


#endif //PAGRAPH_KMERADJNODE_HPP
