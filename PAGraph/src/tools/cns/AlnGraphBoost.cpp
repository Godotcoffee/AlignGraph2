#include <cstdint>
#include <cfloat>
#include <cassert>
#include <string>
#include <queue>
#include <map>
#include <vector>
#include <boost/format.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/utility/value_init.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"

AlnGraphBoost::AlnGraphBoost(const std::string& backbone) {
    // initialize the graph structure with the backbone length + enter/exit
    // vertex
    size_t blen = backbone.length();
    _g = G(blen+2);
    for (size_t i = 0; i < blen+1; i++)
        boost::add_edge(i, i+1, _g);

    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    _enterVtx = *curr++;
    _g[_enterVtx].base = '^';
    _g[_enterVtx].backbone = true;
    for (size_t i = 0; i < blen; i++, ++curr) {
        VtxDesc v = *curr;
        _g[v].backbone = true;
        _g[v].weight = 1;
        _g[v].base = backbone[i];
        _bbMap[v] = v;
    }
    _exitVtx = *curr;
    _g[_exitVtx].base = '$';
    _g[_exitVtx].backbone = true;
}

AlnGraphBoost::AlnGraphBoost(const size_t blen) {
    _g = G(blen+2);
    for (size_t i = 0; i < blen+1; i++)
        boost::add_edge(i, i+1, _g);

    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    _enterVtx = *curr++;
    _g[_enterVtx].base = '^';
    _g[_enterVtx].backbone = true;
    for (size_t i = 0; i < blen; ++i, ++curr) {
        VtxDesc v = *curr;
        _g[v].backbone = true;
        _g[v].weight = 1;
        _g[v].deleted = false;
        _g[v].base = 'N';
        _bbMap[v] = v;
    }
    _exitVtx = *curr;
    _g[_exitVtx].base = '$';
    _g[_exitVtx].backbone = true;
}

void AlnGraphBoost::addAln(dagcon::Alignment& aln, int weight) {
    if (weight <= 0) { return; }

    IndexMap index = boost::get(boost::vertex_index, _g);
    // tracks the position on the backbone
    uint32_t bbPos = aln.start;
    VtxDesc prevVtx = _enterVtx;
    for (size_t i = 0; i < aln.qstr.length(); i++) {
        char queryBase = aln.qstr[i], targetBase = aln.tstr[i];
        assert(queryBase != '.');
        assert(targetBase != '.');
        VtxDesc currVtx = index[bbPos];
        // match
        if (queryBase == targetBase) {
            //_g[_bbMap[currVtx]].coverage++;
            _g[_bbMap[currVtx]].coverage += weight;

            // NOTE: for empty backbones
            _g[_bbMap[currVtx]].base = targetBase;

            //_g[currVtx].weight++;
            _g[currVtx].weight += weight;
            addEdge(prevVtx, currVtx, weight);
            bbPos++;
            prevVtx = currVtx;
        // query deletion
        } else if (queryBase == '-' && targetBase != '-') {
            //_g[_bbMap[currVtx]].coverage++;
            _g[_bbMap[currVtx]].coverage += weight;

            // NOTE: for empty backbones
            _g[_bbMap[currVtx]].base = targetBase;

            bbPos++;
        // query insertion
        } else if (queryBase != '-' && targetBase == '-') {
            // create new node and edge
            VtxDesc newVtx = boost::add_vertex(_g);
            _g[newVtx].base = queryBase;
            //_g[newVtx].weight++;
            _g[newVtx].weight += weight;
            _g[newVtx].backbone = false;
            _g[newVtx].deleted = false;
            _bbMap[newVtx] = bbPos;
            addEdge(prevVtx, newVtx, weight);
            prevVtx = newVtx;
        }
    }
    addEdge(prevVtx, _exitVtx, weight);
}

void AlnGraphBoost::addEdge(VtxDesc u, VtxDesc v, int weight) {
    // Check if edge exists with prev node.  If it does, increment edge counter,
    // otherwise add a new edge.
    InEdgeIter ii, ie;
    bool edgeExists = false;
    for (boost::tie(ii, ie) = boost::in_edges(v, _g); ii != ie; ++ii) {
        EdgeDesc e = *ii;
        if (boost::source(e , _g) == u) {
            // increment edge count
            //_g[e].count++;
            _g[e].count += weight;
            edgeExists = true;
        }
    }
    if (! edgeExists) {
        // add new edge
        std::pair<EdgeDesc, bool> p = boost::add_edge(u, v, _g);
        //_g[p.first].count++;
        _g[p.first].count += weight;
    }
}

void AlnGraphBoost::mergeNodes() {
    std::queue<VtxDesc> seedNodes;
    seedNodes.push(_enterVtx);

    while(true) {
        if (seedNodes.size() == 0)
            break;

        VtxDesc u = seedNodes.front();
        seedNodes.pop();
        mergeInNodes(u);
        mergeOutNodes(u);

        OutEdgeIter oi, oe;
        for (boost::tie(oi, oe) = boost::out_edges(u, _g); oi != oe; ++oi) {
            EdgeDesc e = *oi;
            _g[e].visited = true;
            VtxDesc v = boost::target(e, _g);
            InEdgeIter ii, ie;
            int notVisited = 0;
            for (boost::tie(ii, ie) = boost::in_edges(v, _g); ii != ie; ++ii) {
                if (_g[*ii].visited == false)
                    notVisited++;
            }

            // move onto the boost::target node after we visit all incoming edges for
            // the boost::target node
            if (notVisited == 0)
                seedNodes.push(v);
        }
    }
}

void AlnGraphBoost::mergeInNodes(VtxDesc n) {
    std::map<char, std::vector<VtxDesc>> nodeGroups;
    InEdgeIter ii, ie;
    // Group neighboring nodes by base
    for(boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
        VtxDesc inNode = boost::source(*ii, _g);
        if (out_degree(inNode, _g) == 1) {
            nodeGroups[_g[inNode].base].push_back(inNode);
        }
    }

    // iterate over node groups, merge an accumulate information
    for(auto kvp = nodeGroups.cbegin(); kvp != nodeGroups.end(); ++kvp) {
        std::vector<VtxDesc> nodes = (*kvp).second;
        if (nodes.size() <= 1)
            continue;

        std::vector<VtxDesc>::const_iterator ni = nodes.cbegin();
        VtxDesc an = *ni++;
        OutEdgeIter anoi, anoe;
        boost::tie(anoi, anoe) = boost::out_edges(an, _g);

        // Accumulate out edge information
        for (; ni != nodes.cend(); ++ni) {
            OutEdgeIter oi, oe;
            boost::tie(oi, oe) = boost::out_edges(*ni, _g);
            _g[*anoi].count += _g[*oi].count;
            _g[an].weight += _g[*ni].weight;
        }

        // Accumulate in edge information, merges nodes
        ni = nodes.cbegin();
        ++ni;
        for (; ni != nodes.cend(); ++ni) {
            InEdgeIter ii, ie;
            VtxDesc n = *ni;
            for (boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
                VtxDesc n1 = boost::source(*ii, _g);
                EdgeDesc e;
                bool exists;
                boost::tie(e, exists) = edge(n1, an, _g);
                if (exists) {
                    _g[e].count += _g[*ii].count;
                } else {
                    std::pair<EdgeDesc, bool> p = boost::add_edge(n1, an, _g);
                    _g[p.first].count = _g[*ii].count;
                    _g[p.first].visited = _g[*ii].visited;
                }
            }
            markForReaper(n);
        }
        mergeInNodes(an);
    }
}

void AlnGraphBoost::mergeOutNodes(VtxDesc n) {
    std::map<char, std::vector<VtxDesc>> nodeGroups;
    OutEdgeIter oi, oe;
    for(boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
        VtxDesc outNode = boost::target(*oi, _g);
        if (in_degree(outNode, _g) == 1) {
            nodeGroups[_g[outNode].base].push_back(outNode);
        }
    }

    for(auto kvp = nodeGroups.cbegin(); kvp != nodeGroups.end(); ++kvp) {
        std::vector<VtxDesc> nodes = (*kvp).second;
        if (nodes.size() <= 1)
            continue;

        std::vector<VtxDesc>::const_iterator ni = nodes.cbegin();
        VtxDesc an = *ni++;
        InEdgeIter anii, anie;
        boost::tie(anii, anie) = boost::in_edges(an, _g);

        // Accumulate inner edge information
        for (; ni != nodes.cend(); ++ni) {
            InEdgeIter ii, ie;
            boost::tie(ii, ie) = boost::in_edges(*ni, _g);
            _g[*anii].count += _g[*ii].count;
            _g[an].weight += _g[*ni].weight;
        }

        // Accumulate and merge outer edge information
        ni = nodes.cbegin();
        ++ni;
        for (; ni != nodes.cend(); ++ni) {
            OutEdgeIter oi, oe;
            VtxDesc n = *ni;
            for (boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
                VtxDesc n2 = boost::target(*oi, _g);
                EdgeDesc e;
                bool exists;
                boost::tie(e, exists) = edge(an, n2, _g);
                if (exists) {
                    _g[e].count += _g[*oi].count;
                } else {
                    std::pair<EdgeDesc, bool> p = boost::add_edge(an, n2, _g);
                    _g[p.first].count = _g[*oi].count;
                    _g[p.first].visited = _g[*oi].visited;
                }
            }
            markForReaper(n);
        }
    }
}

void AlnGraphBoost::markForReaper(VtxDesc n) {
    _g[n].deleted = true;
    clear_vertex(n, _g);
    _reaperBag.push_back(n);
}

void AlnGraphBoost::reapNodes() {
    int reapCount = 0;
    std::sort(_reaperBag.begin(), _reaperBag.end());
    std::vector<VtxDesc>::iterator curr = _reaperBag.begin();
    for (; curr != _reaperBag.end(); ++curr) {
        assert(_g[*curr].backbone==false);
        remove_vertex(*curr-reapCount++, _g);
    }
}

const std::string AlnGraphBoost::consensus(int minWeight) {
    // get the best scoring path
    std::vector<AlnNode> path = bestPath();

    // consensus sequence
    std::string cns;

    // track the longest consensus path meeting minimum weight
    int offs = 0, bestOffs = 0, length = 0, idx = 0;
    bool metWeight = false;
    std::vector<AlnNode>::iterator curr = path.begin();
    for (; curr != path.end(); ++curr) {
        AlnNode n = *curr;
        if (n.base == _g[_enterVtx].base || n.base == _g[_exitVtx].base)
            continue;

        cns += n.base;

        // initial beginning of minimum weight section
        if (!metWeight && n.weight >= minWeight) {
            offs = idx;
            metWeight = true;
        } else if (metWeight && n.weight < minWeight) {
        // concluded minimum weight section, update if longest seen so far
            if ((idx - offs) > length) {
                bestOffs = offs;
                length = idx - offs;
            }
            metWeight = false;
        }
        idx++;
    }

    // include end of sequence
    if (metWeight && (idx - offs) > length) {
        bestOffs = offs;
        length = idx - offs;
    }

    return cns.substr(bestOffs, length);
}

void AlnGraphBoost::consensus(std::vector<CnsResult>& seqs, int minWeight, size_t minLen) {
    seqs.clear();

    // get the best scoring path
    std::vector<AlnNode> path = bestPath();

    // consensus sequence
    std::string cns;

    // track the longest consensus path meeting minimum weight
    int offs = 0, idx = 0;
    bool metWeight = false;
    std::vector<AlnNode>::iterator curr = path.begin();
    for (; curr != path.end(); ++curr) {
        AlnNode n = *curr;
        if (n.base == _g[_enterVtx].base || n.base == _g[_exitVtx].base)
            continue;

        cns += n.base;

        // initial beginning of minimum weight section
        if (!metWeight && n.weight >= minWeight) {
            offs = idx;
            metWeight = true;
        } else if (metWeight && n.weight < minWeight) {
        // concluded minimum weight section, add sequence to supplied vector
            metWeight = false;
            CnsResult result;
            result.range[0] = offs;
            result.range[1] = idx;
            size_t length = idx - offs;
            result.seq = cns.substr(offs, length);
            if (length >= minLen) seqs.push_back(result);
        }
        idx++;
    }

    // include end of sequence
    if (metWeight) {
        size_t length = idx - offs;
        CnsResult result;
        result.range[0] = offs;
        result.range[1] = idx;
        result.seq = cns.substr(offs, length);
        if (length >= minLen) seqs.push_back(result);
    }
}

const std::vector<AlnNode> AlnGraphBoost::bestPath() {
    EdgeIter ei, ee;
    for (boost::tie(ei, ee) = edges(_g); ei != ee; ++ei)
        _g[*ei].visited = false;

    std::map<VtxDesc, EdgeDesc> bestNodeScoreEdge;
    std::map<VtxDesc, float> nodeScore;
    std::queue<VtxDesc> seedNodes;

    // start at the end and make our way backwards
    seedNodes.push(_exitVtx);
    nodeScore[_exitVtx] = 0.0f;

    while (true) {
        if (seedNodes.size() == 0)
            break;

        VtxDesc n = seedNodes.front();
        seedNodes.pop();

        bool bestEdgeFound = false;
        float bestScore = -FLT_MAX;
        EdgeDesc bestEdgeD = boost::initialized_value;
        OutEdgeIter oi, oe;
        for(boost::tie(oi, oe) = boost::out_edges(n, _g); oi != oe; ++oi) {
            EdgeDesc outEdgeD = *oi;
            VtxDesc outNodeD = boost::target(outEdgeD, _g);
            AlnNode outNode = _g[outNodeD];
            float newScore, score = nodeScore[outNodeD];
            if (outNode.backbone && outNode.weight == 1) {
                newScore = score - 10.0f;
            } else {
                AlnNode bbNode = _g[_bbMap[outNodeD]];
                newScore = _g[outEdgeD].count - bbNode.coverage*0.5f + score;
            }

            if (newScore > bestScore) {
                bestScore = newScore;
                bestEdgeD = outEdgeD;
                bestEdgeFound = true;
            }
        }

        if (bestEdgeFound) {
            nodeScore[n]= bestScore;
            bestNodeScoreEdge[n] = bestEdgeD;
        }

        InEdgeIter ii, ie;
        for (boost::tie(ii, ie) = boost::in_edges(n, _g); ii != ie; ++ii) {
            EdgeDesc inEdge = *ii;
            _g[inEdge].visited = true;
            VtxDesc inNode = boost::source(inEdge, _g);
            int notVisited = 0;
            OutEdgeIter oi, oe;
            for (boost::tie(oi, oe) = boost::out_edges(inNode, _g); oi != oe; ++oi) {
                if (_g[*oi].visited == false)
                    notVisited++;
            }

            // move onto the target node after we visit all incoming edges for
            // the target node
            if (notVisited == 0)
                seedNodes.push(inNode);
        }
    }

    // construct the final best path
    VtxDesc prev = _enterVtx, next;
    std::vector<AlnNode> bpath;
    while (true) {
        bpath.push_back(_g[prev]);
        if (bestNodeScoreEdge.count(prev) == 0) {
            break;
        } else {
            EdgeDesc bestOutEdge = bestNodeScoreEdge[prev];
            _g[prev].bestOutEdge = bestOutEdge;
            next = boost::target(bestOutEdge, _g);
            _g[next].bestInEdge = bestOutEdge;
            prev = next;
        }
    }

    return bpath;
}

void AlnGraphBoost::printGraph() {
    reapNodes();
    boost::write_graphviz(std::cout, _g,
        make_label_writer(get(&AlnNode::base, _g)),
        make_label_writer(get(&AlnEdge::count, _g)));
}

bool AlnGraphBoost::danglingNodes() {
    VtxIter curr, last;
    boost::tie(curr, last) = boost::vertices(_g);
    bool found = false;
    for (;curr != last; ++curr) {
        if (_g[*curr].deleted)
            continue;
        if (_g[*curr].base == _g[_enterVtx].base || _g[*curr].base == _g[_exitVtx].base)
            continue;

        int indeg = out_degree(*curr, _g);
        int outdeg = in_degree(*curr, _g);
        if (outdeg > 0 && indeg > 0) continue;

        found = true;
    }
    return found;
}

AlnGraphBoost::~AlnGraphBoost(){}
