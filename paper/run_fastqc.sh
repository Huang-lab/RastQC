#!/bin/bash
# Wrapper to run FastQC from compiled source
export PATH="/usr/local/opt/openjdk@11/bin:$PATH"
FASTQC_DIR="/Users/huangk06/Projects/FastQC2/FastQC"
java -Djava.awt.headless=true -Xmx512m \
  -cp "$FASTQC_DIR/bin:$FASTQC_DIR/sam-1.103.jar:$FASTQC_DIR/jbzip2-0.9.jar:$FASTQC_DIR/htsjdk.jar:$FASTQC_DIR/cisd-jhdf5.jar" \
  uk.ac.babraham.FastQC.FastQCApplication "$@"
