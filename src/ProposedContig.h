#ifndef PROPOSED_CONTIG_H
#define PROPOSED_CONTIG_H

#include <deque>
#include <vector>
#include "Contig.h"
#include "KMer.h"

class ProposedContig {
public:
    int length; // Length of the total proposed contig.
    std::deque<Read*>* reads; // Should be reads with relative offsets
    std::deque<int64_t>* readOffsets;
    uint64_t startOffset;
    
    ProposedContig(KMerEdge* initialEdge);
    Contig* exportContig();
    void addKmerForward(KMerEdge* kmerEdge);
    void addKmerBackward(KMerEdge* kmerEdge);
};

#endif
