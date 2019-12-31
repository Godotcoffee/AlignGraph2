#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include "Alignment.hpp"

using namespace dagcon;

///
/// Simple method to reverse complement a sequence.
///
std::string revComp(std::string& seq) {
    const std::string bases = "ACTG";
    std::string::iterator curr = seq.begin();
    for (; curr != seq.end(); ++curr) {
        char& c = *curr;
        c = c == 'T' ? bases[0] :
            c == 'G' ? bases[1] :
            c == 'A' ? bases[2] :
            c == 'C' ? bases[3] : c;
    }
    return std::string(seq.rbegin(), seq.rend());
}

// Set this to false if the alignments are grouped by query.  The parse
// routine will be adjusted to build the alignment graph based on the
// queries.
bool Alignment::groupByTarget = true;

Alignment::Alignment() :
    tlen(0),
    start(0),
    end(0),
    id(""),
    sid(""),
    strand('+'),
    qstr(""),
    tstr("") { }

// Parses blasr m5 output grouped either by target or query.
void parseM5(std::istream& stream, Alignment* aln) {
    std::string line;
    std::getline(stream, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;
    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    // avoids *some* empty lines
    if (fields.size() == 0) return;

    // base query id (without the last '/<coordinates>'), allows us to
    // group properly by query when asked.
    std::string baseQid = fields[0].substr(0,fields[0].find_last_of("/"));
    aln->sid = fields[0];
    aln->id = Alignment::groupByTarget ? fields[5] : baseQid;

    std::istringstream ssLen(Alignment::groupByTarget ? fields[6] : fields[1]);
    ssLen >> aln->tlen;
    std::istringstream ssStart(Alignment::groupByTarget ? fields[7] : fields[2]);
    ssStart >> aln->start;
    aln->start++;

    // the target is always reversed.
    aln->strand = fields[9][0];
    if (aln->strand == '-' && Alignment::groupByTarget) {
        // only need to reverse complement when correcting targets
        aln->qstr = revComp(fields[16]);
        aln->tstr = revComp(fields[18]);
    } else {
        aln->qstr = Alignment::groupByTarget ? fields[16] : fields[18];
        aln->tstr = Alignment::groupByTarget ? fields[18] : fields[16];
    }
}

void parsePre(std::istream& stream, Alignment* aln) {
    std::string line;
    std::getline(stream, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;
    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    // avoids *some* empty lines
    if (fields.size() == 0) return;

    // qid, tid, strand, tlen, tstart, tend, qstr, tstr
    aln->sid = fields[0];
    aln->id = fields[1];
    aln->strand = fields[2][0];

    std::istringstream ssLen(fields[3]);
    ssLen >> aln->tlen;

    std::istringstream ssStart(fields[4]);
    ssStart >> aln->start;

    std::istringstream ssEnd(fields[5]);
    ssEnd >> aln->end;

    aln->qstr = fields[6];
    aln->tstr = fields[7];
}

// default to parsing m5
Alignment::ParseFunc Alignment::parse = parseM5;

std::istream& operator>>(std::istream& instrm, Alignment& data) {
    Alignment::parse(instrm, &data);
    return instrm;
}

std::ostream& operator<<(std::ostream& ostrm, const Alignment& data) {
    ostrm << "target: " << data.id << ", query: " << data.sid;
    ostrm << ", start: " << data.start << ", end: " << data.end;
    ostrm << ", length: " << data.tlen <<  std::endl;
    ostrm << "tstr(50): " << data.tstr.substr(0,50) << std::endl;
    ostrm << "qstr(50): " << data.qstr.substr(0,50) << std::endl;
    return ostrm;
}

Alignment normalizeGaps(Alignment& aln, bool push) {
    // XXX: optimize this
    assert(aln.qstr.length() == aln.tstr.length());
    size_t len = aln.qstr.length();
    std::string qNorm, tNorm;
    qNorm.reserve(len+100);
    tNorm.reserve(len+100);
    std::string qstr = aln.qstr;
    std::string tstr = aln.tstr;

    // convert dots to dashes
    for (size_t i=0; i < len; i++) {
        if ('.' == qstr[i]) qstr[i] = '-';
        if ('.' == tstr[i]) tstr[i] = '-';
    }

    // convert mismatches to indels
    for (size_t i=0; i < len; i++) {
        char qb = qstr[i], tb = tstr[i];
        if (qb != tb && qb != '-' && tb != '-') {
            qNorm += '-';
            qNorm += qb;
            tNorm += tb;
            tNorm += '-';
        } else {
            qNorm += qb;
            tNorm += tb;
        }
    }

    // update length
    assert(qNorm.length() == tNorm.length());
    len = qNorm.length();

    if (push) {
        // push gaps to the right, but not past the end
        for (size_t i=0; i < len-1; i++) {
            // pushing target gaps
            if (tNorm[i] == '-') {
                size_t j = i;
                while (++j < len) {
                    char c = tNorm[j];
                    if (c != '-') {
                        if (c == qNorm[i]) {
                            tNorm[i] = c;
                            tNorm[j] = '-';
                        }
                        break;
                    }
                }
            }

            // pushing query gaps
            if (qNorm[i] == '-') {
                size_t j = i;
                while (++j < len) {
                    char c = qNorm[j];
                    if (c != '-') {
                        if (c == tNorm[i]) {
                            qNorm[i] = c;
                            qNorm[j] = '-';
                        }
                        break;
                    }
                }
            }
        }
    }
    assert(qNorm.length() == tNorm.length());
    assert(len == tNorm.length());

    // generate the final, normalized alignment strings
    Alignment finalNorm;
    finalNorm.id = aln.id;
    finalNorm.sid = aln.sid;
    finalNorm.start = aln.start;
    finalNorm.tlen = aln.tlen;
    finalNorm.strand = aln.strand;
    for (size_t i=0; i < len; i++) {
        if (qNorm[i] != '-' || tNorm[i] != '-') {
            finalNorm.qstr += qNorm[i];
            finalNorm.tstr += tNorm[i];
        }
    }

    return finalNorm;
}

void trimAln(Alignment& aln, int trimLen) {
    int lbases, rbases;
    size_t loffs, roffs;
    auto const len = aln.tstr.length();

    lbases = 0; loffs = 0U;
    while(lbases < trimLen && loffs < len) {
        if (aln.tstr[loffs++] != '-') {
            lbases++;
        }
    }

    rbases = 0;
    roffs = len;
    while (rbases < trimLen && roffs > loffs) {
        if (aln.tstr[--roffs] != '-') {
            rbases++;
        }
    }

    aln.start += lbases;
    aln.qstr = aln.qstr.substr(loffs, roffs - loffs);
    aln.tstr = aln.tstr.substr(loffs, roffs - loffs);
}
