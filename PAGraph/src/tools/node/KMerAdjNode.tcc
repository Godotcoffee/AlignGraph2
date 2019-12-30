//
// Created by amadeus on 19-7-5.
//

template<typename PosType, typename IndexType, typename CountType>
KMerAdjNode<PosType, IndexType, CountType>::KMerAdjNode() {}

template<typename PosType, typename IndexType, typename CountType>
template<typename T>
void KMerAdjNode<PosType, IndexType, CountType>::addWithMutex(std::mutex &mutex, std::vector<T> &store, const T &item) {
    std::lock_guard<std::mutex> lock(mutex);
    store.push_back(item);
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T>
void KMerAdjNode<PosType, IndexType, CountType>::addDirectly(std::vector<T> &store, const T &item) {
    store.push_back(item);
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T>
void KMerAdjNode<PosType, IndexType, CountType>::addItem(std::mutex &mutex, std::vector<T> &store, const T &item,
                                                         bool useMutex) {
    if (useMutex) {
        addWithMutex(mutex, store, item);
    } else {
        addDirectly(store, item);
    }
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T>
void KMerAdjNode<PosType, IndexType, CountType>::addItems(std::mutex &mutex, std::vector<T> &store,
                                                          const std::vector<T> &item, bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(mutex);
        store.insert(std::end(store), std::begin(item), std::end(item));
    } else {
        store.insert(std::end(store), std::begin(item), std::end(item));
    }
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T>
std::size_t KMerAdjNode<PosType, IndexType, CountType>::removeDuplicate(std::vector<T> &array,
                                                                        std::function<bool(const T &, const T &)> cmp,
                                                                        bool shrink) {
    std::size_t reduce = 0;

    if (!array.empty()) {
        std::sort(array.begin(), array.end());

        std::size_t p = 1;
        for (std::size_t k = 1; k < array.size(); ++k) {
            if (!cmp(array[p - 1], array[k])) {
                array[p++] = array[k];
            }
        }

        reduce = array.size() - p;
        array.resize(p);
    }

    if (shrink) {
        array.shrink_to_fit();
    }

    return reduce;
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T1, typename T2>
std::size_t KMerAdjNode<PosType, IndexType, CountType>::cluster(std::vector<T1> &array, std::vector<T2> &counts,
                                                                std::function<bool(const T1 &, const T1 &)> cmp,
                                                                bool shrink) {
    std::size_t reduce = 0;

    std::size_t p = 0;
    auto it1 = array.begin();
    auto it2 = counts.begin();
    for (; it1 != array.end() && it2 != counts.end(); ++it1, ++it2) {
        auto &item = *it1;
        auto &cnt = *it2;

        bool similar = false;

        for (std::size_t i = 0; i < p; ++i) {
            if (cmp(item, array[i])) {
                similar = true;
                counts[i] += cnt;
                break;
            }
        }

        if (!similar) {
            array[p] = item;
            counts[p] = cnt;
            ++p;
        }
    }

    reduce = array.size() - p;
    array.resize(p);

    if (shrink) {
        array.shrink_to_fit();
    }

    return reduce;
}

template<typename PosType, typename IndexType, typename CountType>
template<typename T1, typename T2>
void KMerAdjNode<PosType, IndexType, CountType>::sortWithCount(std::vector<T1> &array, std::vector<T2> &count) {
    std::vector<std::pair<T1, T2>> tmpArr;

    {
        auto it1 = array.begin(), it2 = count.begin();
        for (; it1 != array.end() && it2 != count.end(); ++it1, ++it2) {
            tmpArr.emplace_back(*it1, *it2);
        }
    }

    array.clear();
    count.clear();

    std::sort(tmpArr.begin(), tmpArr.end(), [](const std::pair<T1, T2> &lhs, const std::pair<T1, T2> &rhs) {
        return lhs.first < rhs.first;
    });

    for (auto &p : tmpArr) {
        array.push_back(p.first);
        count.push_back(p.second);
    }
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::reset() {
    _children.clear();
    _positions.clear();
    _posCount.clear();
    _usedFlag.assign(_positions.size(), 0);

    std::vector<IndexType>().swap(_children);
    std::vector<PosType>().swap(_positions);
    std::vector<CountType>().swap(_posCount);
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::addChild(const IndexType &child, bool useMutex) {
    addItem(_mutex, _children, child, useMutex);
}

template<typename PosType, typename IndexType, typename CountType>
std::size_t
KMerAdjNode<PosType, IndexType, CountType>::mergeChild(std::function<bool(const IndexType &, const IndexType &)> cmp,
                                                       bool shrink) {
    return removeDuplicate(_children, cmp, shrink);
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::addPosition(const PosType &pos, bool useMutex) {
    addItem(_mutex, _positions, pos, useMutex);
    addItem(_mutex, _posCount, CountType(1), useMutex);
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::addPosition(const std::vector<PosType> &pos, bool useMutex) {
    addItems(_mutex, _positions, pos, useMutex);
    addItems(_mutex, _posCount, std::vector<CountType>(pos.size(), CountType(1)), useMutex);
}

template<typename PosType, typename IndexType, typename CountType>
std::size_t
KMerAdjNode<PosType, IndexType, CountType>::mergePosition(std::function<bool(const PosType &, const PosType &)> cmp,
                                                          bool shrink) {
    //return removeDuplicate(_positions, cmp, shrink);
    return cluster(_positions, _posCount, cmp, shrink);
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::sortPosition() {
    sortWithCount(_positions, _posCount);
}

template<typename PosType, typename IndexType, typename CountType>
CountType KMerAdjNode<PosType, IndexType, CountType>::getPosCount(std::size_t index, bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        return _posCount[index];
    } else {
        return _posCount[index];
    }
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::resetUsedFlag(bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        _usedFlag.assign(_positions.size(), 0);
    } else {
        _usedFlag.assign(_positions.size(), 0);
    }
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::setFlag(std::size_t index, bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        _usedFlag[index] = true;
    } else {
        _usedFlag[index] = true;
    }
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::unsetFlag(std::size_t index, bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        _usedFlag[index] = false;
    } else {
        _usedFlag[index] = false;
    }
}

template<typename PosType, typename IndexType, typename CountType>
bool KMerAdjNode<PosType, IndexType, CountType>::getFlag(std::size_t index, bool useMutex) const {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        return _usedFlag[index];
    } else {
        return _usedFlag[index];
    }
}

template<typename PosType, typename IndexType, typename CountType>
const std::vector<IndexType> &KMerAdjNode<PosType, IndexType, CountType>::getAllChild() const {
    return _children;
}

template<typename PosType, typename IndexType, typename CountType>
const std::vector<PosType> &KMerAdjNode<PosType, IndexType, CountType>::getAllPositions() const {
    return _positions;
}

template<typename PosType, typename IndexType, typename CountType>
const std::vector<CountType> &KMerAdjNode<PosType, IndexType, CountType>::getAllCount() const {
    return _posCount;
}

template<typename PosType, typename IndexType, typename CountType>
void KMerAdjNode<PosType, IndexType, CountType>::resetAll(bool useMutex) {
    if (useMutex) {
        std::lock_guard<std::mutex> lock(_mutex);
        reset();
    } else {
        reset();
    }
}

