//
// Created by amadeus on 19-12-2.
//

template<typename Selector, typename Filter, typename T>
void PAlgorithm::filterHelper(std::vector<T> &contents, Selector s, Filter f) {
    std::size_t p = 0;
    for (std::size_t i = 0; i < contents.size(); ++i) {
        if (f(s(contents[i]))) {
            contents[p++] = contents[i];
        }
    }
    contents.resize(p);
}

template<typename Filter>
void PAlgorithm::filterSuccessors(const PABruijnGraph::PANode &paNode,
                                  std::vector<std::pair<PABruijnGraph::PANode, int>> &successors, Filter f) {
    /*std::size_t p = 0;
    for (std::size_t i = 0; i < successors.size(); ++i) {
        auto &successor = successors[i];
        if (f(paNode, successor)) {
            successors[p++] = successors[i];
        }
    }
    successors.resize(p);*/

    filterHelper(
            successors,
            [](const std::pair<PABruijnGraph::PANode, int> &i) -> const std::pair<PABruijnGraph::PANode, int> &
            { return i; },
            [&](const std::pair<PABruijnGraph::PANode, int> &i) { return f(paNode, i); });
}

template<typename Filter>
void
PAlgorithm::classifySuccessors(std::vector<std::pair<PABruijnGraph::PANode, int>> &results,
                               const PABruijnGraph::PANode &paNode, PABruijnGraph::PosType deviation, double errorRate,
                               std::pair<std::int64_t, std::int64_t> ctgRange, bool canLeap, double leapMin,
                               const PositionMapper &ctgMapper, Filter filter) {
    auto rootPos = paNode.getPosition();

    std::vector<std::pair<PABruijnGraph::PANode, int>> originSuccessors;

    _pPagraph->successors(originSuccessors, paNode, deviation, errorRate);

    filterSuccessors(paNode, originSuccessors, filter);

    std::vector<std::size_t> amazing;    // best: Ctg and Ref
    std::vector<std::size_t> excellent;  // better: go to ctg or go to ref with ctg
    std::vector<std::size_t> great;      // good: only ref or only ctg
    std::vector<std::size_t> skip;       // ctg skip: ctg jump


    for (std::size_t i = 0; i < originSuccessors.size(); ++i) {
        auto &successor = originSuccessors[i];
        auto pos = successor.first.getPosition();
        auto check = PABruijnGraph::checkPosition(rootPos, pos, successor.second, deviation, errorRate);

        auto leap = pos.first != 0
                    && (pos.first < ctgRange.first || pos.first >= ctgRange.second);

        if (leap) {
            auto dual = ctgMapper.singleToDual(pos.first);
            if (dual.second > ctgMapper.size(dual.first) * leapMin) {
                continue;
            }
        }

        if (!canLeap && leap) {
            continue;
        }

        if (check == PABruijnGraph::MatchGrade::Amazing || leap) {
            amazing.emplace_back(i);
        } else if (check == PABruijnGraph::MatchGrade::Excellent) {
            excellent.emplace_back(i);
        } else if (check == PABruijnGraph::MatchGrade::Good) {
            great.emplace_back(i);
        } else if (canLeap && check == PABruijnGraph::MatchGrade::Skip) {
            skip.emplace_back(i);
        }
    }

    auto &chosenOne = !amazing.empty() ? amazing : (!excellent.empty() ? excellent : (!great.empty() ? great : skip));

    for (auto idx : chosenOne) {
        results.push_back(originSuccessors[idx]);
    }
}


template<typename ParentFilter>
PAlgorithm::NodeStatus PAlgorithm::walkStraight(const std::pair<PABruijnGraph::PANode, int> &startNode,
                                                std::vector<std::pair<PABruijnGraph::PANode, int>> &path, int deviation,
                                                double errorRate, std::pair<std::int64_t, std::int64_t> ctgRange,
                                                std::size_t hasSize, std::size_t splitSize, double splitMin,
                                                const PositionMapper &ctgMapper, ParentFilter parentFilter, size_t limitation) {

    std::set<PABruijnGraph::PANode::UniqueType> uniqueTable;
    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> ctgPosTable;
    resetCtgPosTable(ctgPosTable);

    auto &paNode = startNode.first;
    auto dist = startNode.second;

    std::size_t nowSize = dist;

    path.emplace_back(paNode, dist);

    if (paNode.getPosition().first != 0 &&
        (paNode.getPosition().first < ctgRange.first || paNode.getPosition().first >= ctgRange.second)) {
        return NodeStatus::Leap;
    }

    insertCtgPosHelper(ctgPosTable, path.front().first.getPosition().first);

    uniqueTable.insert(paNode.getUniqueHelper());

    auto filter = [&](const PABruijnGraph::PANode &paNode,
                      const std::pair<PABruijnGraph::PANode, int> &successor) -> bool {

        return parentFilter(paNode, successor)
               &&
               uniqueTable.count(successor.first.getUniqueHelper()) == 0
               &&
               (successor.first.getPosition().first == 0 ||
                PABruijnGraph::isEdgeSimilar(
                        paNode.getPosition(),
                        successor.first.getPosition(),
                        successor.second,
                        deviation, errorRate).first ||
                !existCtgPos(ctgPosTable, successor.first.getPosition().first));
    };

    std::vector<std::pair<PABruijnGraph::PANode, int>> successors;

    for (;;) {
        successors.clear();
        classifySuccessors(successors, path.back().first, deviation, errorRate, ctgRange,
                                             (hasSize + nowSize) >= splitSize, splitMin, ctgMapper, filter);

        if (successors.empty()) {
            return NodeStatus::End;
        }

        if (successors.size() > 1) {
            return NodeStatus::Branch;
        }

        auto &front = successors.front();
        auto checkItem = front.first.getUniqueHelper();

        uniqueTable.insert(checkItem);
        insertCtgPosHelper(ctgPosTable, front.first.getPosition().first);

        path.push_back(front);
        nowSize += front.second;

        auto lastCtgPos = path.back().first.getPosition().first;

        if (lastCtgPos != 0 && (lastCtgPos < ctgRange.first || lastCtgPos >= ctgRange.second)) {
            return NodeStatus::Leap;
        }

        if (limitation > 0 && path.size() >= limitation) {
            return NodeStatus::Limit;
        }
    }
}

template<typename ParentFilter>
PAlgorithm::TravelSequence PAlgorithm::graphTravel(const PABruijnGraph::PANode &startNode, int deviation,
                                                   double errorRate,
                                                   std::pair<std::int64_t, std::int64_t> ctgRange, std::size_t hasSize,
                                                   std::size_t splitSize, double splitMin,
                                                   const PositionMapper &ctgMapper, ParentFilter parentFilter) {

    std::pair<PABruijnGraph::PosType, PABruijnGraph::PosType> ctgTravelPosTable;
    resetCtgPosTable(ctgTravelPosTable);

    std::set<PABruijnGraph::PANode::UniqueType> travelUniqueTable;

    std::vector<std::pair<PABruijnGraph::PANode, int>> seq;
    std::size_t nowSize = _pPagraph->getKmerSize();

    std::vector<std::pair<PABruijnGraph::PANode, int>> path;
    std::vector<decltype(path)> paths;

    std::pair<PABruijnGraph::PANode, int> chosenOne = {startNode, static_cast<int>(_pPagraph->getKmerSize())};
    insertCtgPosHelper(ctgTravelPosTable, startNode.getPosition().first);

    auto filter = [&](const PABruijnGraph::PANode &paNode,
                      const std::pair<PABruijnGraph::PANode, int> &successor) -> bool {

        return parentFilter(paNode, successor)
               &&
               travelUniqueTable.count(successor.first.getUniqueHelper()) == 0
               &&
               (successor.first.getPosition().first == 0 ||
                PABruijnGraph::isEdgeSimilar(
                        paNode.getPosition(),
                        successor.first.getPosition(),
                        successor.second,
                        deviation, errorRate).first ||
                !existCtgPos(ctgTravelPosTable, successor.first.getPosition().first));
    };

    walkStraight(chosenOne, path, deviation, errorRate, ctgRange, hasSize + nowSize, splitSize, splitMin, ctgMapper, filter);
    paths.emplace_back(path);

    std::size_t chosenOndIdx = 0;

    std::vector<std::pair<PABruijnGraph::PANode, int>> successors;

    for (;;) {
        //path.clear();
        //std::cout << "Walk straight " << cnt << std::endl;

        auto &chosenPath = paths[chosenOndIdx];

        //std::cout << "Walk straight Done" << cnt++ << std::endl;

        for (auto &p : chosenPath) {
            seq.push_back(p);
            travelUniqueTable.insert(p.first.getUniqueHelper());
            nowSize += p.second;
        }

        insertCtgPosHelper(ctgTravelPosTable, chosenPath);

        auto &lastNode = seq.back().first;

        auto lastCtgPos = lastNode.getPosition().first;

        if (lastCtgPos != 0 && (lastCtgPos < ctgRange.first || lastCtgPos >= ctgRange.second)) {
            break;
        }

        successors.clear();
        classifySuccessors(successors, lastNode, deviation, errorRate, ctgRange,
                                             (hasSize + nowSize) >= splitSize, splitMin, ctgMapper, filter);

        std::vector<std::pair<std::size_t, std::size_t>> leap;
        std::vector<std::pair<std::size_t, std::size_t>> branch;
        std::vector<std::pair<std::size_t, std::size_t>> tips;

        paths.clear();

        for (std::size_t i = 0; i < successors.size(); ++i) {
            auto &successor = successors[i];

            path.clear();

            auto status = walkStraight(successor, path, deviation, errorRate, ctgRange, hasSize + nowSize, splitSize, splitMin, ctgMapper, filter);
            paths.emplace_back(path);

            if (status == NodeStatus::Leap) {
                leap.emplace_back(i, path.size());
            } else if (status == NodeStatus::End) {
                tips.emplace_back(i, path.size());
            } else {
                branch.emplace_back(i, path.size());
            }
        }

        if (leap.empty() && tips.empty() && branch.empty()) {
            break;
        }

        if (!leap.empty()) {
            chosenOndIdx = leap.front().first;
            //chosenOne = successors[chosenOndIdx];
        } else if (!branch.empty()) {
            std::size_t chosen = 0;
            for (std::size_t i = 1; i < branch.size(); ++i) {
                auto &paNode1 = successors[branch[i].first].first;
                auto &paNode2 = successors[branch[chosen].first].first;
                if (paNode1.getAbundance() > paNode2.getAbundance()) {
                    chosen = i;
                }
            }
            chosenOndIdx = branch[chosen].first;
            //chosenOne = successors[chosenOndIdx];
        } else {
            std::size_t chosen = 0;
            for (std::size_t i = 1; i < tips.size(); ++i) {
                if (tips[i].second > tips[chosen].second) {
                    chosen = i;
                }
            }
            chosenOndIdx = tips[chosen].first;
            //chosenOne = successors[chosenOndIdx];
        }
    }

    return seq;
}

template<typename Functor>
void PAlgorithm::searchPANode(std::vector<std::pair<PABruijnGraph::ANode, size_t>> &aNodes,
                              std::vector<PABruijnGraph::PANode> &result, bool onlyFirst, Functor f) {

    std::set<PABruijnGraph::PANode::UniqueType> unique;

    for (auto &node : aNodes) {
        auto &aNode = node.first;
        for (decltype(aNode.size()) i = 0; i < aNode.size(); ++i) {
            auto pos = aNode[i];
            auto dualPos1 = _pCtgMapper->singleToDual(static_cast<std::size_t>(pos.first));
            auto dualPos2 = _pRefMapper->singleToDual(static_cast<std::size_t>(pos.second));

            auto paNode = PABruijnGraph::PANode(aNode, i);
            auto uniqueHelper = paNode.getUniqueHelper();

            if (!aNode.isUsed(i) && unique.count(uniqueHelper) == 0 &&
                f(node.second, dualPos1.first, dualPos1.second, dualPos2.first, dualPos2.second)) {

                result.push_back(paNode);
                unique.insert(uniqueHelper);
            }
        }
        if (!result.empty() && onlyFirst) {
            break;
        }
    }
}
