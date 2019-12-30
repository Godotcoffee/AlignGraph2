//
// Created by amadeus on 19-12-26.
//

#ifndef PAGRAPH_FILEKMERITERATOR_HPP
#define PAGRAPH_FILEKMERITERATOR_HPP


#include "IKmerIterator.hpp"

class FileKmerIterator : public IKmerIterator {
private:
    std::string _filePath;
    std::size_t _k;
public:
    explicit FileKmerIterator(std::string path);

    size_t kSize() const override { return _k; };

    void iterate(std::function<void(std::uint64_t)> f, unsigned threadNum) const override;
};


#endif //PAGRAPH_FILEKMERITERATOR_HPP
