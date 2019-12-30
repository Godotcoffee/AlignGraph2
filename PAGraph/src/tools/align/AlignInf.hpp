#include <utility>

//
// Created by amadeus on 19-5-5.
//

#ifndef PAGRAPH_ALIGNINF_HPP
#define PAGRAPH_ALIGNINF_HPP

#include <cstddef>
#include <string>
#include <vector>

/**
 * Structure of align information,
 * with a default comparision using score (Larger).
 */
class AlignInf {
private:
    std::string _queryName;
    std::string _refName;
    std::size_t _score;
    std::size_t _queryBegin;
    std::size_t _queryEnd;
    std::size_t _refBegin;
    std::size_t _refEnd;
    bool _forward;

    std::vector<bool> _queryDiff;
    std::vector<bool> _refDiff;

public:
    explicit AlignInf();

    explicit AlignInf(
            std::string queryName,
            std::string refName,
            std::size_t score,
            std::size_t queryBegin,
            std::size_t queryEnd,
            std::size_t refBegin,
            std::size_t refEnd,
            bool forward,
            std::vector<bool> queryDiff,
            std::vector<bool> refDiff);

    bool operator<(const AlignInf &rhs) const;

    const std::string &getQueryName() const;

    void setQueryName(const std::string &queryName);

    const std::string &getRefName() const;

    void setRefName(const std::string &refName);

    std::size_t getScore() const;

    void setScore(std::size_t score);

    std::size_t getQueryBegin() const;

    void setQueryBegin(std::size_t queryBegin);

    std::size_t getQueryEnd() const;

    void setQueryEnd(std::size_t queryEnd);

    std::size_t getRefBegin() const;

    void setRefBegin(std::size_t refBegin);

    std::size_t getRefEnd() const;

    void setRefEnd(std::size_t refEnd);

    bool isForward() const;

    void setForward(bool forward);

    const std::vector<bool> &getQueryDiff() const;

    const std::vector<bool> &getRefDiff() const;
};

#endif //PAGRAPH_ALIGNINF_HPP
