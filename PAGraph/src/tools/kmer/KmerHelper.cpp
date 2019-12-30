//
// Created by amadeus on 19-12-30.
//

#include "KmerHelper.hpp"

void KmerHelper::kmer2Code(std::vector<std::uint64_t> &codes, const std::string &seq, std::size_t kmerSize) {
    auto len = (std::size_t) seq.size();
    std::uint64_t code = 0;

    std::uint64_t indexMask = (1UL << (kmerSize * 2)) - 1;

    for (std::size_t i = 0; i < kmerSize && i < len; ++i) {
        code = (code << 2U) | acgt(seq[i]);
    }
    if (len >= kmerSize) {
        codes.emplace_back(code);
    }

    for (std::size_t i = kmerSize; i < len; ++i) {
        code = ((code << 2U) | acgt(seq[i])) & indexMask;

        codes.emplace_back(code);
    }
}

std::string KmerHelper::code2Kmer(std::uint64_t code, std::size_t kmerSize) {
    std::string kmer((std::size_t) kmerSize, 'a');
    const char *table = "ACGT";
    auto idx = code;
    for (decltype(kmerSize) i = 0; i < kmerSize; ++i) {
        kmer[kmerSize - 1 - i] = table[idx & 0x3U];
        idx >>= 2U;
    }

    return kmer;
}
