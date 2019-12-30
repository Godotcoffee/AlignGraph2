//
// Created by amadeus on 19-5-5.
//

#include <iostream>
#include "CompressedSeq.hpp"

CompressedSeq::CompressedSeq(const std::string &seq) {
    _length = seq.size();
    if (_length == 0) {
        return;
    }
    unsigned buffer = 0;
    for (decltype(_length) i = 0; i < _length; ++i) {
        unsigned encode = 0;
        switch (seq[i]) {
            case 'A': case 'a': default:
                encode = 0; break;
            case 'C': case 'c':
                encode = 1; break;
            case 'G': case 'g':
                encode = 2; break;
            case 'T': case 't':
                encode = 3; break;
        }
        buffer = buffer | (encode << ((i & 0x3U) * 2U));
        // i % 4 == 3
        if ((i & 0x3U) == 0x3U) {
            _data.push_back((unsigned char) (buffer & 0xFFU));
            buffer = 0;
        }
    }

    // length % 4 != 0
    if (_length & 0x3U) {
        _data.push_back((unsigned char) (buffer & 0xFFU));
    }
}

void CompressedSeq::partDecode(std::string &buf, std::size_t idx, bool forward) const {
    const char *table = forward ? "ACGT" : "TGCA";

    buf.clear();
    unsigned code = _data[idx];
    
    for (int i = 0; i < 4; ++i) {
        buf.push_back(table[code & 0x3U]);
        code >>= 2U;
    }
}

std::size_t CompressedSeq::size() const {
    return _length;
}

std::string CompressedSeq::toString(bool forward) const {
    const char *table = forward ? "ACGT" : "TGCA";
    std::string seq;
    seq.resize(_length);
    std::size_t index = 0;
    if (_length > 0) {
        unsigned decode = 0;
        for (std::size_t i = 0; i < _length; ++i) {
            // i % 4 == 0
            if ((i & 0x3U) == 0) {
                decode = _data[index++];
            }
            //seq.push_back(table[decode & 0x3]);
            seq[forward ? i : _length - 1 - i] = table[decode & 0x3U];
            decode >>= 2U;
        }
    }
    return seq;
}

char CompressedSeq::baseAt(std::size_t idx, bool forward) const {
    if (idx >= _length) {
        return 'N';
    }
    auto trueIdx = forward ? idx / 4 : (_length - 1 - idx) / 4;
    auto offset = forward ? idx % 4 : (_length - 1 - idx) % 4;

    std::string buf;
    partDecode(buf, trueIdx, forward);

    return buf[offset];
}
