#include <stdint.h>
#include <stdlib.h>
#include <deque>
#include "Contig.h"
#include "KMer.h"
#include "Read.h"
#include "ProposedContig.h"


// 
//   ###########[SRC MATCH]########
//              [EXT MATCH]###################
//

ProposedContig::ProposedContig(KMerEdge* initialEdge) {
    this->reads.push_back(initialEdge->src->source);
    this->length = initialEdge->src->source->length; // Start as the length of one read
    this->startOffset = 0;
    
    this->readOffsets.push_back(initialEdge->src->offset - initialEdge->ext->offset); // The next read should be overlayed ahead by the src match less the ext (although in current implementation it is always 0) 
    this->reads.push_back(initialEdge->ext->source);
    this->length += initialEdge->src->offset - initialEdge->ext->offset + initialEdge->ext->source->length - initialEdge->src->source->length; // Add the overhang length, account for reads of different length
}

void ProposedContig::addKmerForward(KMerEdge* kmerEdge) {
    this->readOffsets.push_back(kmerEdge->src->offset - kmerEdge->ext->offset); // The next read should be overlayed ahead by the src match less the ext (although in current implementation it is always 0) 
    this->reads.push_back(kmerEdge->ext->source);
    this->length += kmerEdge->src->offset - kmerEdge->ext->offset + kmerEdge->ext->source->length - kmerEdge->src->source->length; // Add the overhang length, account for reads of different length
}

void ProposedContig::addKmerBackward(KMerEdge* kmerEdge) {
    /*
    this->reads.push_front(kmerEdge->src->source);
    this->readOffsets.push_front(kmerEdge->distance);
    this->startOffset += kmerEdge->distance;
    this->length += kmerEdge->distance;*/
}

Contig* ProposedContig::exportContig() {
    Read* pop = this->reads.front();
    int64_t popOffset;
    uint64_t accumulatedOffset = 0;
    char* sequence = (char*) malloc(this->length + 1);
    while (!this->reads.empty()) {
        this->reads.pop_front();
        for (int i = 0; i < (pop->length); i++) {
            sequence[i + accumulatedOffset] = pop->sequence[i];
            //sequence[i + accumulatedOffset] = 'a';
        }
        pop = this->reads.front();
        popOffset = this->readOffsets.front();
        this->readOffsets.pop_front();
        accumulatedOffset = accumulatedOffset + popOffset;
        this->swapReads.push_front(pop);
        this->swapReadOffsets.push_front(popOffset);
    }
    while (!this->swapReads.empty()) {
        pop = this->swapReads.front();
        this->swapReads.pop_front();
        this->reads.push_front(pop);
    }
    while (!this->swapReadOffsets.empty()) {
        popOffset = this->swapReadOffsets.front();
        this->swapReadOffsets.pop_front();
        this->readOffsets.push_front(popOffset);
    }
    return new Contig(sequence);
}
