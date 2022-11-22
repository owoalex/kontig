#include <stdint.h>
#include <stdlib.h>
#include <deque>
#include "Contig.h"
#include "KMer.h"
#include "ProposedContig.h"

ProposedContig::ProposedContig(KMerEdge* initialEdge) {
    this->reads.push_front(initialEdge->src->source);
    this->readOffsets.push_front(0);
    this->reads.push_front(initialEdge->ext->source);
    this->readOffsets.push_front(initialEdge->distance);
    //this->length = initialEdge->src->length + initialEdge->distance;
    this->endOffset = initialEdge->distance;
    this->startOffset = 0;
}

void ProposedContig::addKmerForward(KMerEdge* kmerEdge) {
    this->reads.push_front(kmerEdge->ext->source);
    this->readOffsets.push_front(kmerEdge->distance);
    this->endOffset += kmerEdge->distance;
    this->length += kmerEdge->distance;
}

void ProposedContig::addKmerBackward(KMerEdge* kmerEdge) {
    this->reads.push_back(kmerEdge->src->source);
    this->readOffsets.push_back(kmerEdge->distance);
    this->startOffset += kmerEdge->distance;
    this->length += kmerEdge->distance;
}

Contig* ProposedContig::exportContig() {
    Read* pop;
    int64_t popOffset;
    char* sequence = (char*) malloc(this->length);
    while (!this->reads.empty()) {
        pop = this->reads.front();
        popOffset = this->readOffsets.front();
        this->reads.pop_front();
        this->readOffsets.pop_front();
        this->swapReads.push_front(pop);
        this->swapReadOffsets.push_front(popOffset);
    }
    while (!this->swapReads.empty()) {
        pop = this->swapReads.front();
        popOffset = this->swapReadOffsets.front();
        this->swapReads.pop_front();
        this->swapReadOffsets.pop_front();
        this->reads.push_front(pop);
        this->readOffsets.push_front(popOffset);
    }
    return new Contig(sequence);
}
