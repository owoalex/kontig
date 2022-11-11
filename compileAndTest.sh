#!/bin/bash
if make; then
    build/kontig -i sample_17_R1_001.fastq
else
    echo "BUILD FAILED"
fi
