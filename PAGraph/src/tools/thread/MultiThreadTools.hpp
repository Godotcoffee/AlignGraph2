//
// Created by amadeus on 19-5-10.
//

#ifndef PAGRAPH_MULTITHREADTOOLS_HPP
#define PAGRAPH_MULTITHREADTOOLS_HPP

#include <thread>
#include <vector>
#include <memory>
#include "IMultiThread.hpp"

/**
 * Multiple threads helper
 */
class MultiThreadTools {
private:
    template<typename T, typename Functor>
    static void traversalHelper(T &obj, unsigned threadNum, Functor &fun);

    template<typename T, typename Functor>
    static void traversalHelper(const T &obj, unsigned threadNum, Functor &fun);
public:
    template<typename T, typename Functor>
    static void multiTraversal(IMultiThread<T> &obj, unsigned threadNum, Functor fun);

    template<typename T, typename Functor>
    static void multiTraversal(const IMultiThread<T> &obj, unsigned threadNum, Functor fun);

    template<typename T, typename Functor>
    static void multiTraversal(std::vector<T> &obj, unsigned threadNum, Functor fun);

    template<typename T, typename Functor>
    static void multiTraversal(const std::vector<T> &obj, unsigned threadNum, Functor fun);
};

#include "MultiThreadTools.tcc"

#endif //PAGRAPH_MULTITHREADTOOLS_HPP
