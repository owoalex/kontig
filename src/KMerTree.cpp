#include <stdint.h>
#include <stdlib.h>
#include "KMer.h"
#include "KMerTree.h"

KMerTree::KMerTree() {
    this->root = new KMerSortingTreeNode();
    this->root->high = nullptr;
    this->root->low = nullptr;
    this->root->values = nullptr;
}

void KMerTree::insertKmer(KMer* kmer) {
    KMerSortingTreeNode* cnode = this->root;
    for (int j = 0; j < 64; j++) {
        if ((kmer->quickRef >> j) & 0b1) {
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
    KMerLinkedList* ll = cnode->values;
    if (ll == nullptr) {
        ll = new KMerLinkedList();
        ll->next = nullptr;
        ll->value = kmer;
        cnode->values = ll;
    } else {
        while (ll->next != nullptr) {
            ll = ll->next;
        }
        ll->next = new KMerLinkedList();
        ll->next->next = nullptr;
        ll->next->value = kmer;
    }
}
