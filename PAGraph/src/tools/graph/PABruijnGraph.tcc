//
// Created by amadeus on 19-7-10.
//

template<typename T>
void PABruijnGraph::sampleSequence(std::vector<std::pair<PABruijnGraph::KmerIndex, std::size_t>> &result,
                                   const std::string &seq, std::size_t outerSample, T functor) const {
    auto len = seq.size();
    auto index = buildCode(seq);
    std::int64_t lastOne = -1;

    if (len >= _kmerSize) {
        for (std::size_t i = 0; i < len - _kmerSize + 1; ++i) {
            auto idx1 = index[i];

            if (functor(static_cast<std::size_t>(i))) {
                auto check = searchDenseIndex(idx1);
                if (check.second) {
                    if (lastOne < 0 || i - static_cast<std::size_t>(lastOne) >= outerSample) {
                        result.emplace_back(check.first, i);
                        lastOne = i;
                    }
                }
            }
        }
    }
}