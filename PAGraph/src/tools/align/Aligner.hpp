#include <utility>

//
// Created by alice on 18-11-7.
//

#ifndef PAGRAPH_ALIGNER_HPP
#define PAGRAPH_ALIGNER_HPP


#include <cstddef>
#include <unordered_map>
#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <memory>
#include "seq/ISeqDatabase.hpp"
#include "seq/SeqInf.hpp"
#include "AlignReference.hpp"
#include "AlignInf.hpp"
#include "IAlignDatabase.hpp"
#include "thread/MultiThreadTools.hpp"
#include "ParseAlignTools.hpp"

class Aligner {
public:
    typedef std::int64_t pos_t;
    typedef std::pair<std::int64_t, pos_t> ref_pos_t;
private:
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReadDB;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pContigDB;
    std::shared_ptr<const ISeqDatabase<SeqInf>> _pReferenceDB;

    std::shared_ptr<const IAlignDatabase<AlignInf>> _pReadToCtgDB;
    std::shared_ptr<const IAlignDatabase<AlignInf>> _pReadToRefDB;
    std::shared_ptr<const IAlignDatabase<AlignInf>> _pCtgToRefDB;

    // Arguments
    int     _readToCtgTopK = 3;
    int     _readToRefTopK = 3;
    int     _ctgToRefTopK = -1;
    double  _readToCtgRatio = 0.5;
    double  _readToRefRatio = 0.4;
    double  _ctgToRefRatio = 0.1;
    double  _ctgToRefTotalRatio = 0.3;
    std::size_t
            _covFilter = 1;
    std::size_t _ctgToRefMinLen = 2500;

    /**
     * Store a more index for convenience
     */
    class ExAlignInf {
    private:
        AlignInf _baseAlignInf;
        std::size_t _refIdx;        // Store an index for reference
    public:
        explicit ExAlignInf(AlignInf alignInf, std::size_t refIdx) : _baseAlignInf(std::move(alignInf)), _refIdx(refIdx) {}

        const AlignInf &getBaseAlignInf() {
            return _baseAlignInf;
        }

        std::size_t getRefIdx() {
            return _refIdx;
        }

        bool operator<(const ExAlignInf &rhs) const {
            return _baseAlignInf < rhs._baseAlignInf;
        }
    };

    /* Align */
    std::vector<std::vector<ExAlignInf>> _readToContig;
    std::vector<std::vector<ExAlignInf>> _readToRef;
    std::vector<std::vector<ExAlignInf>> _contigToRef;

    /* Coverage */
    std::vector<std::vector<std::size_t>> _refCov;

    /* Alignments between contigs and references */
    AlignReference _alignReference;

    /* Flags for references */
    std::vector<bool> _refFilterFlag;

    /* Flags for references */
    std::vector<bool> _ctgFilterFlag;
    std::vector<bool> _ctgFilterForward;

    bool _isInit;

    void init(bool debug = false);

    void mergeAlignInf();

    static void
    mergeAlignInfHelper(std::vector<std::vector<ExAlignInf>> &align, const IAlignDatabase<AlignInf> &alignDB,
                        const ISeqDatabase<SeqInf> &queryDB, const ISeqDatabase<SeqInf> &refDB);

    static void
    covInfHelper(std::vector<std::vector<std::size_t>> &align,
                               const IAlignDatabase<AlignInf> &alignDB,
                               const ISeqDatabase<SeqInf> &queryDB, const ISeqDatabase<SeqInf> &refDB);

    bool queryContig(std::size_t contigIdx, std::size_t contigPos, bool forward,
                     std::vector<std::pair<ref_pos_t, ref_pos_t>> &positions) const;

    template <typename T>
    void parseToCtg(T functor, unsigned threadNum = 8);

    template <typename T>
    void parseToRef(T functor, unsigned threadNum = 8);

    static void flipPosition(std::size_t &left, std::size_t &right, std::size_t length);

    static void clearHelper(std::vector<std::vector<std::pair<ref_pos_t, ref_pos_t>>> &vec) {
        for (auto &v : vec) {
            v.clear();
        }
    };
public:
    explicit Aligner(
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReadDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pContigDB,
            std::shared_ptr<const ISeqDatabase<SeqInf>> pReferenceDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToContigDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pReadToRefDB,
            std::shared_ptr<const IAlignDatabase<AlignInf>> pContigToRefDB
    );

    /**
     * Simple align from contig to reference.
     * Flags can be set by calling setCtgFilter() or setRefFilter().
     * @return
     */
    void simpleAlign();

    /**
     * Add 0 to unalign position of contigs.
     * It can be controlled by setUseCtgAlign(). Contig with False flag will be ignored.
     */
    void addExtraPosition();

    /**
     * Change the flags of contig.
     * @param ctgIdx
     * @param forward
     * @param use
     * @return
     */
    bool setUseCtgAlign(const std::string &ctgName, bool forward, bool use);

    template<typename T>
    void runReadToCtg(T functor,
                      unsigned threadNum = 8,
                      bool debug = false);

    template<typename T>
    void runReadToRef(T functor,
                      unsigned threadNum = 8,
                      bool debug = false);

    void setReadToCtgTopK(int readToCtgTopK);

    void setReadToRefTopK(int readToRefTopK);

    void setCtgToRefTopK(int ctgToRefTopK);

    void setReadToCtgRatio(double readToCtgRatio);

    void setReadToRefRatio(double readToRefRatio);

    void setCtgToRefRatio(double ctgToRefRatio);

    void setCtgToRefTotalRatio(double ctgToRefTotalRatio);

    void setCtgToRefMinLen(std::size_t ctgToRefMinLen);

    void setCovFilter(std::size_t cov);

    bool setRefFilter(const std::string &refName, bool accepted);

    void clearRefFilter(bool accepted = true);

    bool setCtgFilter(const std::string &refName, bool forward, bool accepted);

    void clearCtgFilter(bool accepted = true);
};

#include "Aligner.tcc"


#endif //PAGRAPH_ALIGNER_HPP
