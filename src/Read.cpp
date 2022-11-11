#include <stdint.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#include "Read.h"
#include "KMer.h"

Read::Read(char* name, char* sequence, char* qualities) {
    this->length = 0;
    while (sequence[this->length] != '\0') {
        this->length++;
    }
    this->sequence = sequence;
    this->name = name;
    this->qualities = (uint8_t*) malloc(this->length);
    for (int i = 0; i < this->length; i++) {
        this->qualities[i] = qualities[i] - 0x21;
    }
}

KMerSet* Read::generateKmers(int kmerLength) {
    this->kmers = (KMerSet*) malloc(sizeof(KMerSet));
    this->kmers->length = 0;
    if (length > kmerLength) {
        this->kmers->length = (length - kmerLength);
        this->kmers->kmers = (KMer**) malloc(sizeof(KMer*) * this->kmers->length);
        for (int i = 0; i < this->kmers->length; i++) {
            this->kmers->kmers[i] = new KMer(this->sequence + i, this->qualities + i, this, i, kmerLength);
        }
    }
    return this->kmers;
}
