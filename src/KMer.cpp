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
    this->quickRef = 0;
    for (int i = 0; i < this->length && i < 32; i++) {
        this->quickRef = this->quickRef << 2;
        switch (this->sequence[i]) {
            case 'A':
                this->quickRef = this->quickRef | 0b00;
                break;
            case 'T':
                this->quickRef = this->quickRef | 0b01;
                break;
            case 'G':
                this->quickRef = this->quickRef | 0b10;
                break;
            case 'C':
                this->quickRef = this->quickRef | 0b11;
                break;
        }
    }
    //this->quickRef = ;
}

bool KMer::isEqual(KMer* other) {
    if (other->length != this->length) {
        return false; // We don't currently match mismatched kmers at all
    }
    for (int i = 0; i < this->length; i++) {
        if (this->sequence[i] != other->sequence[i]) {
            return false;
        }
    }
    return true;
}

uint64_t KMer::getMatchQuality(KMer* other) {
    uint64_t q = 0;
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
