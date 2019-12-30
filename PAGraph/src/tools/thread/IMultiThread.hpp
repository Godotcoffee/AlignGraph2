//
// Created by amadeus on 19-5-10.
//

#ifndef PAGRAPH_IMULTITHREAD_HPP
#define PAGRAPH_IMULTITHREAD_HPP

template<typename T>
class IMultiThread {
public:
    /**
     * @return Total size of data.
     */
    virtual std::size_t size() const = 0;

    /**
     * const getter
     * @param id
     * @return
     */
    virtual const T &operator[](std::size_t id) const = 0;

    /**
     * getter
     * @param id
     * @return
     */
    virtual T &operator[](std::size_t id) = 0;
};

#endif //PAGRAPH_IMULTITHREAD_HPP
