# Identifying-Intron-Junctions-from-Aligned-RNA-Seq-Reads
Script for detecting and counting splice junctions from RNA-Seq alignment data (SAM files). Developed as coursework for MSc Bioinformatics 2024/25 at the University of Glasgow.



This script identifies splice junctions from RNA-Seq alignments stored in a SAM file and reports those that fall within annotated gene regions.

**Functionality:**

* Parses RNA-Seq alignments from a SAM file, ignoring header lines.
* Detects split reads based on the presence of the `N` operator in the CIGAR string.
* Extracts genomic coordinates for each splice junction.
* Ensures only uniquely aligned reads (`NH:i:1`) are considered.
* Counts the number of reads supporting each unique junction.
* Matches junctions to gene locations provided in a separate annotation file.
* Outputs a summary listing:

  * Gene ID
  * Junction start and end positions
  * Number of supporting reads
* Inserts a blank line between results for each gene for clarity.

**Usage:**

```
python3 script.py alignments.sam annotation.txt
```

**Output:**

* A `.txt` file listing all junctions found within annotated gene regions.
* Each line corresponds to one junction.
* Blank lines separate gene sections.

**Requirements:**

* Python 3
* No external libraries needed; uses only built-in Python functions.

