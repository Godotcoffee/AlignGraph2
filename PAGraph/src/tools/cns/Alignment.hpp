#pragma once

#include <stdint.h>
#include <iostream>

///
/// Super-simple alignment representation.  Represents an alignment between two
/// PacBio reads, one of which we're trying to correct.  The read to correct
/// may be either the target or the query, depending on how the alignment was
/// done.
///
namespace dagcon {
class Alignment {
public:
    typedef void (*ParseFunc)(std::istream&, Alignment* aln);

    // May correct the target or the query, default is target
    static bool groupByTarget;

    // length of the sequence we are trying to correct
    uint32_t tlen;

    // conforming offsets are 1-based
    uint32_t start;

    uint32_t end;

    // ID of the read we're trying to correct (target)
    std::string id;

    // ID of the supporting read (query)
    std::string sid;

    char strand;

    // query and target strings must be equal length
    std::string qstr;
    std::string tstr;

    Alignment();

    static ParseFunc parse;
};
}

std::istream& operator>>(std::istream& instrm, dagcon::Alignment& data);
std::ostream& operator<<(std::ostream& ostrm, const dagcon::Alignment& data);

void parseM5(std::istream& stream, dagcon::Alignment* aln);

void parsePre(std::istream& stream, dagcon::Alignment* aln);

/// Simplifies the alignment by normalizing gaps.  Converts mismatches into
/// indels ...
///      query: CAC        query:  C-AC
///             | |  --->          |  |
///     target: CGC       target:  CG-C
///
/// Shifts equivalent gaps to the right in the reference ...
///      query: CAACAT        query: CAACAT
///             | | ||  --->         |||  |
///     target: C-A-AT       target: CAA--T
///
/// Shifts equivalent gaps to the right in the read ...
///      query: -C--CGT       query: CCG--T
///              |  | |  --->        |||  |
///     target: CCGAC-T      target: CCGACT
/// Allow optional gap pushing, some aligners may not need it and I'd like
/// to get rid of it anyway.
dagcon::Alignment normalizeGaps(dagcon::Alignment& aln, bool push=true);

void trimAln(dagcon::Alignment& aln, int trimLen=50);

std::string revComp(std::string& seq);
