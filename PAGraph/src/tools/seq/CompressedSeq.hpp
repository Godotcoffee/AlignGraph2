//
// Created by amadeus on 19-5-5.
//

#ifndef PAGRAPH_COMPRESSEDSEQ_HPP
#define PAGRAPH_COMPRESSEDSEQ_HPP

#include <vector>
#include <string>


/**
 * Store the sequence with 2 bits per base.
 * Note: It can't store the all types of base except ACGT, which will be transfer to 'A'.
 */
class CompressedSeq {
private:
    std::vector<unsigned char> _data;
    std::size_t _length;

    void partDecode(std::string &buf, std::size_t idx, bool forward = true) const;
public:
    explicit CompressedSeq(const std::string &seq);

    std::size_t size() const;

    std::string toString(bool forward = true) const;

    char baseAt(std::size_t, bool forward = true) const;
};


#endif //PAGRAPH_COMPRESSEDSEQ_HPP
