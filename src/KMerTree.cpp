#include <stdint.h>
#include <stdlib.h>
#include "KMer.h"
#include "KMerTree.h"
#include "KMerTreeBuilder.h"

int KMerTree::powi (int base, unsigned int exp){
    int res = 1;
    while (exp) {
        if (exp & 1)
            res *= base;
        exp >>= 1;
        base *= base;
    }
    return res;
}

KMerTree::KMerTree(KMerTreeBuilder* from) {
    size_t tree_size = sizeof(void*) * (powi(2,(16+1))-1);
    uint64_t kmer_arrays = 0;
    size_t data_size = sizeof(void*) * kmer_arrays;
    heap = malloc(tree_size);
}

KMerSet* KMerTree::getMatchingKmers() {
    //sizeof(void*)
}
