//
// Created by amadeus on 19-12-30.
//

#ifndef PAGRAPH_KMERHELPER_HPP
#define PAGRAPH_KMERHELPER_HPP

#include <vector>
#include <cstdint>
#include <string>

class KmerHelper {
private:
    static unsigned acgt(char ch) {
        int n = 0;
        switch (ch) {
            case 'A': case 'a': default:
                n = 0;
                break;
            case 'C': case 'c':
                n = 1;
                break;
            case 'G': case 'g':
                n = 2;
                break;
            case 'T': case 't':
                n = 3;
                break;
        }
        return n;
    }
public:
    static void kmer2Code(std::vector<std::uint64_t> &codes, const std::string &seq, std::size_t kmerSize);

    static std::string code2Kmer(std::uint64_t code, std::size_t kmerSize);
};


#endif //PAGRAPH_KMERHELPER_HPP
