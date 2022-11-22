#ifndef KMER_H
#define KMER_H

#include <vector>

class Read;

struct KMerEdge;

class KMer {
public:
    int length;
    char* sequence;
    uint8_t* qualities;
    Read* source;
    int offset;
    uint64_t quickref;
    uint64_t basegc; // Denotes if base is G/C
    uint64_t basedir; // Denotes if base is A/G
    std::vector<KMerEdge*> kmerEdgesForward;
    std::vector<KMerEdge*> kmerEdgesBackward;
    
    KMer(char* sequence, uint8_t* qualities, Read* source, int offset, int length);
    bool isEqual(KMer* other);
    uint16_t getMismatchedBases(KMer* other);
    int64_t getMatchQuality(KMer* other);
};

struct KMerSet {
    int length;
    KMer** kmers;
};

struct KMerEdge {
    KMer* src; // Source
    KMer* ext; // Extension
    int64_t weight;
};

#endif
