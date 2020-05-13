//
// Created by amadeus on 19-5-13.
//

#ifndef PAGRAPH_PASSEMBLY_HPP
#define PAGRAPH_PASSEMBLY_HPP

#include "PABruijnGraph.hpp"
#include "PAlgorithm.hpp"
#include "position/PositionMapper.hpp"

class PAssembly {
private:
    template<typename T>
    static void combatSeq(
            const std::vector<PAlgorithm::TravelSequence> &seqs,
            const PositionMapper &mapper,
            std::size_t start,
            bool forward,
            T functor);

public:
    static std::set<std::pair<std::string, bool>> testTravel5(
            const std::string &outDir,
            const std::string &prefix,
            std::shared_ptr<PABruijnGraph> pGraph,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pRefDB,
            std::shared_ptr<const PositionMapper> pCtgMapper,
            std::shared_ptr<const PositionMapper> pRefMapper,
            const std::set<std::pair<std::string, bool>> &ctgSet,
            std::size_t deviation,
            double errorRate,
            double startSplit,
            std::size_t minLen,
            unsigned threadNum);
};

#include "PAssembly.tcc"

#endif //PAGRAPH_PASSEMBLY_HPP
