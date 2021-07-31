# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | r1 |
| 1294_S1_L008_R2_001.fastq.gz | r2 |
| 1294_S1_L008_R3_001.fastq.gz | r3 |
| 1294_S1_L008_R4_001.fastq.gz | r4 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.

r1 histogram: /home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-first/r1_Mean_Qscores.png
r2 histogram: /home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-first/r2_Mean_Qscores.png
r3 histogram: /home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-first/r3_Mean_Qscores.png
r4 histogram: /home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-first/r4_Mean_Qscores.png

    2. ```A Qscore of 20 is my cutoff. Reading through the literature this seems to be the standard. Signifies a 1% chance of an incorrect base call.```
    3. ```r2: 3976613 Ns
          r3: 3328051 Ns```
        ```zcat 1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l```
        ```zcat 1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep 'N' | wc -l```
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
