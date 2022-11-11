#include <stdint.h>
#include <stdlib.h>
#include "Read.h"

class Contig {
public:
    int length;
    char* sequence;
    uint16_t* qualities;
    
    Contig(char* sequence, uint16_t* qualities) {
        this->length = 0;
        while (sequence[this->length] != '\0') {
            this->length++;
        }
        this->sequence = sequence;
        this->qualities = qualities;
    }
};
