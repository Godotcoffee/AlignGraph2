//
// Created by amadeus on 19-7-5.
//

template<typename T, typename Functor>
void MultiThreadTools::traversalHelper(T &obj, unsigned threadNum, Functor &fun) {
    std::vector<std::thread> threads;

    for (decltype(threadNum) t = 0; t < threadNum; ++t) {
        threads.emplace_back([&, t] {
            auto totalSeq = obj.size();
            for (auto i = static_cast<decltype(totalSeq)>(t); i < totalSeq; i += threadNum) {
                fun(t, i, obj[i]);
            }
        });
    }

    for (auto &thread : threads) {
        thread.join();
    }
}

template<typename T, typename Functor>
void MultiThreadTools::traversalHelper(const T &obj, unsigned threadNum, Functor &fun) {
    std::vector<std::thread> threads;

    for (decltype(threadNum) t = 0; t < threadNum; ++t) {
        threads.emplace_back([&, t] {
            auto totalSeq = obj.size();
            for (auto i = static_cast<decltype(totalSeq)>(t); i < totalSeq; i += threadNum) {
                fun(t, i, obj[i]);
            }
        });
    }

    for (auto &thread : threads) {
        thread.join();
    }
}

template<typename T, typename Functor>
void MultiThreadTools::multiTraversal(IMultiThread<T> &obj, unsigned threadNum, Functor fun) {
    traversalHelper(obj, threadNum, fun);
}

template<typename T, typename Functor>
void MultiThreadTools::multiTraversal(const IMultiThread<T> &obj, unsigned threadNum, Functor fun) {
    traversalHelper(obj, threadNum, fun);
}

template<typename T, typename Functor>
void MultiThreadTools::multiTraversal(std::vector<T> &obj, unsigned threadNum, Functor fun) {
    traversalHelper(obj, threadNum, fun);
}

template<typename T, typename Functor>
void MultiThreadTools::multiTraversal(const std::vector<T> &obj, unsigned threadNum, Functor fun) {
    traversalHelper(obj, threadNum, fun);
}
