# Live Virus Analysis

​	This folder contains the code we used to analyze NGS sequencing of SARS-CoV-2 3-chymotrypsin-like protease (3CL <sub>pro</sub>) to count how frequently each amino acid appears at each residue. The NGS sequencing was done by tiling multiple reads across the entire 3CL<sub>pro</sub>, which we take advantage of by filtering reads based on the presence of a known primer that identifies which segment that read belongs to.  

### Prerequisites

Python 3

[Biopython](https://biopython.org/)

[Pandas](https://pandas.pydata.org/)

[EMBOSS](http://emboss.sourceforge.net/)

### Set up

​	Running the pipeline requires the following elements:

- /config folder:
  - conditions.csv : FASTQ files must be named in accordance with the condition names in this file.
  - segment_primers.csv : contains the primers that identify which segment of 3CL <sub>pro</sub> the read belongs to.
-  /reads folder:
  - 1 FASTQ file for each condition.
  - ref3CL.fasta : contains the nucleotide sequence that codes for the reference  3CL <sub>pro</sub>.

### Running the pipeline

​	The pipeline is run with the command `python Live_virus_analysis.py` in the scripts folder.

### Pipeline overview

1. The reads are searched for the segment identifying primers and categorized based on which segment they cover. Identical reads are grouped together to speed alignment.
2. Categorized and grouped reads are saved to a fasta file for each segment, containing the read and the number of times that read was seen.
3. The EMBOSS water pairwise sequence alignment tool is used to align the reads in each segment to the reference nucleotide sequence.
4. Reads containing insertions and deletions are counted and discarded.
5. Reads containing ambiguous residues are counted and discarded.
6. Aligned reads are placed in frame and translated to produce amino acid sequences.
7. Amino acid sequences are analyzed to produce a count of how many times each amino acid appeared at each residue.

### Output

​	A CSV for each condition, containing the number of times each amino acid was seen at each residue.

​	A CSV containing summary statistics, breaking down the number of reads in total, as well as the number of reads that contained indels or X's. 