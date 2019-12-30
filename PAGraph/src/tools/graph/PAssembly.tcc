
template<typename T>
void PAssembly::combatSeq(const std::vector<PAlgorithm::TravelSequence> &seqs, const PositionMapper &mapper,
                          std::size_t start, bool forward, T functor) {
    std::size_t nextIdx = start * 2 + (forward ? 0 : 1);
    std::size_t nextPos = 0;

    std::set<std::size_t> ctgSet;
    ctgSet.emplace(nextIdx);

    while (true) {
        bool go = functor(nextIdx / 2, nextIdx % 2 == 0, nextPos);
        if (!go) {
            break;
        }

        if (seqs[nextIdx].empty() || seqs[nextIdx].back().first.getPosition().first == 0) {
            break;
        }

        auto lastCtgPos = mapper.singleToDual(seqs[nextIdx].back().first.getPosition().first);
        nextPos = lastCtgPos.second;
        nextIdx = (std::abs(lastCtgPos.first) - 1) * 2 + (lastCtgPos.first > 0 ? 0 : 1);

        if (ctgSet.count(nextIdx) > 0) {
            break;
        }

        ctgSet.emplace(nextIdx);
    }
}
