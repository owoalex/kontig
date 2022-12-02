#include <stdint.h>
#include <stdlib.h>
#include "KMer.h"
#include "KMerTreeBuilder.h"

KMerTreeBuilder::KMerTreeBuilder() {
    this->root = new KMerSortingTreeNode();
    this->root->high = nullptr;
    this->root->low = nullptr;
    this->root->values = nullptr;
}

void KMerTreeBuilder::insertKmer(KMer* kmer) {
    KMerSortingTreeNode* cnode = this->root;
    for (int j = 0; j < 64; j++) {
        if ((kmer->quickref >> j) & 0b1) {
            if (cnode->high == nullptr) {
                cnode->high = new KMerSortingTreeNode();
                cnode->high->high = nullptr;
                cnode->high->low = nullptr;
                cnode->high->values = nullptr;
            }
            cnode = cnode->high;
        } else {
            if (cnode->low == nullptr) {
                cnode->low = new KMerSortingTreeNode();
                cnode->low->high = nullptr;
                cnode->low->low = nullptr;
                cnode->low->values = nullptr;
            }
            cnode = cnode->low;
        }
    }
    
    if (cnode->values == nullptr) {
        cnode->values = new KMerLinkedList();
        cnode->values->next = nullptr;
        cnode->values->value = kmer;
    } else {
        KMerLinkedList* ll = new KMerLinkedList();
        ll->next = cnode->values;
        ll->value = kmer;
        cnode->values = ll;
    }
}
