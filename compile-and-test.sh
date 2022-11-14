#!/bin/bash
if make; then
    build/kontig genmap -i sample-h3n2.fastq --ksize 64 -o sample-h3n2.kmf
    #build/kontig genmap --input sample-staph.fastq --ksize 32 -o sample-staph.kmf
else
    echo "BUILD FAILED"
fi
