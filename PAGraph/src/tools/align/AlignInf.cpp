//
// Created by amadeus on 19-7-5.
//

#include "AlignInf.hpp"

AlignInf::AlignInf() :
        _queryName(""),
        _refName(""),
        _score(0),
        _queryBegin(0),
        _queryEnd(0),
        _refBegin(0),
        _refEnd(0),
        _forward(false) {}

AlignInf::AlignInf(std::string queryName, std::string refName, std::size_t score, std::size_t queryBegin,
                   std::size_t queryEnd, std::size_t refBegin, std::size_t refEnd, bool forward,
                   std::vector<bool> queryDiff, std::vector<bool> refDiff) :
        _queryName(std::move(queryName)),
        _refName(std::move(refName)),
        _score(score),
        _queryBegin(queryBegin),
        _queryEnd(queryEnd),
        _refBegin(refBegin),
        _refEnd(refEnd),
        _forward(forward),
        _queryDiff(std::move(queryDiff)),
        _refDiff(std::move(refDiff)) {}

bool AlignInf::operator<(const AlignInf &rhs) const {
    return _score > rhs._score;
}

const std::string &AlignInf::getQueryName() const {
    return _queryName;
}

void AlignInf::setQueryName(const std::string &queryName) {
    AlignInf::_queryName = queryName;
}

const std::string &AlignInf::getRefName() const {
    return _refName;
}

void AlignInf::setRefName(const std::string &refName) {
    AlignInf::_refName = refName;
}

std::size_t AlignInf::getScore() const {
    return _score;
}

void AlignInf::setScore(std::size_t score) {
    AlignInf::_score = score;
}

std::size_t AlignInf::getQueryBegin() const {
    return _queryBegin;
}

void AlignInf::setQueryBegin(std::size_t queryBegin) {
    AlignInf::_queryBegin = queryBegin;
}

std::size_t AlignInf::getQueryEnd() const {
    return _queryEnd;
}

void AlignInf::setQueryEnd(std::size_t queryEnd) {
    AlignInf::_queryEnd = queryEnd;
}

std::size_t AlignInf::getRefBegin() const {
    return _refBegin;
}

void AlignInf::setRefBegin(std::size_t refBegin) {
    AlignInf::_refBegin = refBegin;
}

std::size_t AlignInf::getRefEnd() const {
    return _refEnd;
}

void AlignInf::setRefEnd(std::size_t refEnd) {
    AlignInf::_refEnd = refEnd;
}

bool AlignInf::isForward() const {
    return _forward;
}

void AlignInf::setForward(bool forward) {
    AlignInf::_forward = forward;
}

const std::vector<bool> &AlignInf::getQueryDiff() const {
    return _queryDiff;
}

const std::vector<bool> &AlignInf::getRefDiff() const {
    return _refDiff;
}
