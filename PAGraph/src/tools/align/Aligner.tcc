//
// Created by amadeus on 19-7-5.
//

template<typename T>
void Aligner::runReadToCtg(T functor, unsigned int threadNum, bool debug){
    //init(debug);
    //parsePosition(functor, refIdx, 1, 1, 1);
    if (debug) { std::cout << "start parsing read to contig" << std::endl; }
    parseToCtg(functor, threadNum);
    //parsePosition(functor, refIdx, 0, 0, 0);
}

template<typename T>
void Aligner::runReadToRef(T functor, unsigned int threadNum, bool debug){
    //init(debug);
    //parsePosition(functor, refIdx, 1, 1, 1);
    if (debug) { std::cout << "start parsing read to contig" << std::endl; }
    parseToRef(functor, threadNum);
    //parsePosition(functor, refIdx, 0, 0, 0);
}

template<typename T>
void Aligner::parseToCtg(T functor, unsigned int threadNum) {
    auto &readDB = *_pReadDB;
    auto &contigDB = *_pContigDB;

    MultiThreadTools::multiTraversal
            (readDB, threadNum, [&](decltype(threadNum) t, decltype(readDB.size()) i, const SeqInf &seq) {
                auto &toContig = _readToContig[i];
                std::size_t readLen = seq.size();

                bool useful[2] = {false, false};

                std::vector<std::vector<std::pair<ref_pos_t, ref_pos_t>>> fwdPositions(readLen);
                std::vector<std::vector<std::pair<ref_pos_t, ref_pos_t>>> revPositions(readLen);

                /* Reads to contigs */
                decltype(_readToCtgTopK) ctgTopKCnt = 0;
                for (auto &exAlign : toContig) {
                    if (_readToCtgTopK >= 0 && ctgTopKCnt >= _readToCtgTopK) { break; }

                    auto &align = exAlign.getBaseAlignInf();
                    auto contigIdx = exAlign.getRefIdx();

                    if (!_ctgFilterFlag[contigIdx]) { continue; }

                    auto readBegin = align.getQueryBegin();
                    auto readEnd = align.getQueryEnd();

                    if ((readEnd - readBegin) * 1.0 / readLen < _readToCtgRatio) { continue; }

                    auto isForward = align.isForward();

                    auto contigBegin = align.getRefBegin();
                    auto contigEnd = align.getRefEnd();

                    auto &contigInf = contigDB[contigIdx];
                    auto contigLen = contigInf.size();

                    // Out of range
                    if (contigEnd >= contigLen || contigBegin >= contigLen) {
                        continue;
                    }

                    if (!isForward) {
                        flipPosition(readBegin, readEnd, readLen);
                    }

                    auto &queryDiff = align.getQueryDiff();
                    auto &refDiff = align.getRefDiff();

                    for (int ii = 0; ii < 2; ++ii) {
                        if ((ii == 0) == _ctgFilterForward[contigIdx]) {
                            auto &positions = isForward ? fwdPositions : revPositions;

                        //if (_alignReference.isAlign(contigIdx, ii == 0) || true) {

                            ParseAlignTools::exactAlign(readBegin, contigBegin, ii == 0, queryDiff, refDiff,
                                                        [&](std::size_t curRead, std::size_t curContig) {
                                                            if (curRead >= readBegin && curRead < positions.size()) {
                                                                auto result = queryContig(contigIdx,
                                                                                          curContig, ii == 0,
                                                                                          positions[curRead]);
                                                                useful[isForward ? 0 : 1] =
                                                                        useful[isForward ? 0 : 1] || result;
                                                            }
                                                        });


                        //}
                        }
                        isForward = !isForward;
                        flipPosition(readBegin, readEnd, readLen);
                        flipPosition(contigBegin, contigEnd, contigLen);
                    }

                    ++ctgTopKCnt;
                }

                functor(i, fwdPositions, revPositions, useful[0], useful[1], t);
            });
}

template<typename T>
void Aligner::parseToRef(T functor, unsigned int threadNum) {
    auto &readDB = *_pReadDB;

    MultiThreadTools::multiTraversal

            (readDB, threadNum, [&](decltype(threadNum) t, decltype(readDB.size()) i, const SeqInf &seq) {

                std::size_t readLen = seq.size();

                std::vector<std::vector<std::pair<ref_pos_t, ref_pos_t>>> fwdPositions(readLen);
                std::vector<std::vector<std::pair<ref_pos_t, ref_pos_t>>> revPositions(readLen);

                auto &toRef = _readToRef[i];
                bool useful[2] = {false, false};

                decltype(_readToRefTopK) refTopKCnt = 0;
                for (auto &exAlign : toRef) {
                    if (_readToRefTopK >= 0 && refTopKCnt >= _readToRefTopK) { break; }

                    auto refIdx = exAlign.getRefIdx();

                    if (!_refFilterFlag[refIdx]) { continue; }

                    auto &align = exAlign.getBaseAlignInf();
                    auto readBegin = align.getQueryBegin();
                    auto readEnd = align.getQueryEnd();

                    if ((readEnd - readBegin) * 1.0 / readLen < _readToRefRatio) { continue; }

                    auto isForward = align.isForward();
                    auto refBegin = align.getRefBegin();
                    auto &queryDiff = align.getQueryDiff();
                    auto &refDiff = align.getRefDiff();

                    std::size_t maxCov = 0;

                    for (std::size_t pos = align.getRefBegin(); pos < align.getRefEnd(); ++pos) {
                        if (pos >= _refCov[refIdx].size()) { break; }
                        maxCov = std::max(maxCov, _refCov[refIdx][pos]);
                    }

                    if (maxCov < _covFilter) { continue; }

                    if (!isForward) {
                        flipPosition(readBegin, readEnd, readLen);
                    }

                    auto &positions = isForward ? fwdPositions : revPositions;

                    useful[isForward ? 0 : 1] = true;

                    ParseAlignTools::exactAlign(readBegin, refBegin, true, queryDiff, refDiff, [&](std::size_t readCur, std::size_t refCur) {
                        if (readCur >= readBegin && readCur < positions.size()) {
                            positions[readCur].emplace_back(
                                    ref_pos_t(0, 0),
                                    ref_pos_t(refIdx + 1, refCur));
                        }

                    });

                    ++refTopKCnt;
                }

                functor(i, fwdPositions, revPositions, useful[0], useful[1], t);
            });
}
