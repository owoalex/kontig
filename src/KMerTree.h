#ifndef KMER_TREE_H
#define KMER_TREE_H

class KMer;

struct KMerLinkedList {
    KMer* value;
    KMerLinkedList* next;
};

struct KMerSortingTreeNode {
    KMerSortingTreeNode* high;
    KMerSortingTreeNode* low;
    KMerLinkedList* values;
};

class KMerTree {
public:
    KMerSortingTreeNode* root;
    
    KMerTree();
    void insertKmer(KMer* kmer);
};

#endif
