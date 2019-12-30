//
// Created by amadeus on 19-7-5.
//

template<typename T>
void
exactAlign2(std::size_t queryBegin, std::size_t refBegin, bool forward, const std::vector<int> &diff,
                            std::size_t diffSize, T functor) {
    // forward 表示ref是否反向互补
    auto curRef = refBegin;
    auto curQuery = queryBegin;
    auto curDiff = forward ? diff.begin() : diff.end();
    for (decltype(diffSize) jj = 0; jj < diffSize; ++jj) {
        if ((forward && curDiff == diff.end()) || (!forward && curDiff == diff.begin())) {
            functor(curQuery, curRef);
            ++curRef;
            ++curQuery;
        } else {
            auto &diffItem = forward ? *curDiff : *std::prev(curDiff);
            auto val = diffItem >= 0 ? diffItem : -diffItem - 1;
            if (jj < (forward ? val : diffSize - 1 - val)) {
                functor(curQuery, curRef);
                ++curRef;
                ++curQuery;
            } else {
                if (diffItem >= 0) {
                    ++curRef;
                } else {
                    functor(curQuery, curRef);
                    ++curQuery;
                }

                if (forward) {
                    ++curDiff;
                } else {
                    --curDiff;
                }
            }
        }
    }
}


template<typename T>
void
ParseAlignTools::exactAlign(std::size_t queryBegin, std::size_t refBegin, bool forward,
                            const std::vector<bool> &queryDiff, const std::vector<bool> &refDiff, T functor) {
    if (queryDiff.empty()) {
        return;
    }
    // forward 表示ref是否反向互补
    auto curRef = refBegin;
    auto curQuery = queryBegin;

    for (decltype(queryDiff.size()) jj = 0; jj < queryDiff.size(); ++jj) {
        auto idx = forward ? jj : queryDiff.size() - jj - 1;
        if (!(queryDiff[idx] ^ refDiff[idx])) {
            functor(curQuery, curRef);
            ++curRef;
            ++curQuery;
        } else {
            if (queryDiff[idx]) {
                ++curRef;
            } else {
                functor(curQuery, curRef);
                ++curQuery;
            }
        }
    }
}
