#include <stdint.h>
#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#include "Contig.h"

Contig::Contig(char* sequence) {
    this->length = 0;
    while (sequence[this->length] != '\0') {
        this->length++;
    }
    this->sequence = sequence;
}
