#ifndef KMER_TREE_BUILDER_H
#define KMER_TREE_BUILDER_H

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

class KMerTreeBuilder {
public:
    KMerSortingTreeNode* root;
    
    KMerTreeBuilder();
    void insertKmer(KMer* kmer);
};

#endif
