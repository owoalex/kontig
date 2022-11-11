#ifndef READ_H
#define READ_H

struct KMerSet;

class Read {
public:
    int length;
    char* sequence;
    char* name;
    uint8_t* qualities;
    KMerSet* kmers;
    
    Read(char* name, char* sequence, char* qualities);
    KMerSet* generateKmers(int kmerLength);
};

#endif
