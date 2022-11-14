#!/bin/bash
if make; then
    build/kontig genmap --fastq sample_17_R1_001.fastq --map sample.kmf
else
    echo "BUILD FAILED"
fi
