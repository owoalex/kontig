#ifndef READ_H
#define READ_H

#include <vector>

struct KMerSet;
struct KMerEdge;

class Read {
public:
    int length;
    char* sequence;
    char* name;
    uint8_t* qualities;
    KMerSet* kmers;
    std::vector<KMerEdge*> kmerEdgesForward;
    std::vector<KMerEdge*> kmerEdgesBackward;
    uint64_t usedNTimes;
    
    Read(char* name, char* sequence, char* qualities);
    KMerSet* generateKmers(int kmerLength);
};

#endif
