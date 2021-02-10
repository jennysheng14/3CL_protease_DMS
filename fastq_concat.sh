#!/bin/bash
# Given fastq.gz files demultiplexed by Basespace, this script will concatenate
# all the files into a new folder and run UMI-tools on the reverse reads
# to pick out the UMI and protein barcodes
# as well as on the forward reads to pick out well-barcode 1 and well-barcode
# 2. Downstream processing will be done via python.

# -input -i refers to the directory where all the Basespace fastq.gz files are.
# -output -o refers to the directory where the concatenated files should reside.
# -file -f refers to the .txt file that has a list of all the indices of the
# files in the experiment.

while getopts i:o:f: option
do
case "${option}"
in
i) INPUT=${OPTARG};;
o) OUTPUT=${OPTARG};;
f) FILE=${OPTARG};;
esac
done

mkdir -p $OUTPUT
# while read file; do cp -r $INPUT/"$file"_*/"$file"*R1* $OUTPUT; done < $FILE
while read file; do cat $INPUT/"$file"_*/"$file"*R1* > \
$OUTPUT/"$file"_R1.fastq.gz; done < $FILE
