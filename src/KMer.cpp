#include <stdint.h>
#include <stdlib.h>
#include "KMer.h"
#include "Read.h"

KMer::KMer(char* sequence, uint8_t* qualities, Read* source, int offset, int length) {
    this->length = length;
    this->sequence = sequence;
    this->qualities = qualities;
    this->source = source;
    this->offset = offset;
    this->basegc = 0;
    this->quickref = 0;
    this->basedir = 0;
    for (int i = 0; i < this->length && i < 64; i++) {
        this->basegc = this->basegc << 1;
        this->basedir = this->basedir << 1;
        switch (this->sequence[i]) {
            case 'A':
                this->basegc = this->basegc | 0b0;
                this->basedir = this->basedir | 0b1;
                break;
            case 'T':
                this->basegc = this->basegc | 0b0;
                this->basedir = this->basedir | 0b0;
                break;
            case 'G':
                this->basegc = this->basegc | 0b1;
                this->basedir = this->basedir | 0b1;
                break;
            case 'C':
                this->basegc = this->basegc | 0b1;
                this->basedir = this->basedir | 0b0;
                break;
        }
    }
    this->quickref = this->basegc ^ this->basedir;
}

bool KMer::isEqual(KMer* other) {
    if (other->length != this->length) {
        return false; // We don't currently match mismatched kmers at all
    }
    if (other->basegc != this->basegc) {
        return false;
    }
    if (other->basedir != this->basedir) {
        return false;
    }
    if (this->length > 64) { // We only need to do sequence match on long kmers
        for (int i = 0; i < this->length; i++) {
            if (this->sequence[i] != other->sequence[i]) {
                return false;
            }
        }
    }
    return true;
}

uint16_t KMer::getMismatchedBases(KMer* other) {
    uint64_t gcdiff = other->basegc ^ this->basegc;
    uint64_t dirdiff = other->basedir ^ this->basedir;
    uint64_t diff = gcdiff | dirdiff;
    return __builtin_popcountll(diff);
    
}

int64_t KMer::getMatchQuality(KMer* other) {
    int64_t q = 0;
    if (other->length != this->length) {
        return q; // We don't currently match mismatched kmers at all
    }
    for (int i = 0; i < this->length; i++) {
        if (this->sequence[i] == other->sequence[i]) {
            q += (this->qualities[i] * other->qualities[i]);
        } else {
            q -= (this->qualities[i] * other->qualities[i]);
        }
    }
    return q;
}
