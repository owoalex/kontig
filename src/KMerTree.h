#ifndef KMER_TREE_H
#define KMER_TREE_H

class KMer;
class KMerTreeBuilder;
struct KMerSet;

class KMerTree {
private:
    int powi (int base, unsigned int exp);
public:
    void* heap;
    
    KMerTree(KMerTreeBuilder* from);
    KMerSet* getMatchingKmers();
};

#endif
