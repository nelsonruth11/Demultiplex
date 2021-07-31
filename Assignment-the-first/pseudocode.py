problem: data from 24 different dual-matched libraries, some indexes have hopped and don't match the proper read

solution: loop thru files one line at a time, find indexes that don't match, send all 4 files to the hopped bin
    -output: 48 FASTQ files: (index 1 and 2) + read 1 OR read 2
              2 FASTQ files: hopped read-pairs of index 1 OR index 2
              2 FASTQ files: undermined/low quality of index 1 OR index 2

open R1, R2, R3, R4:

def demultiplexer:
    """Demultiplexes FASTQ files, accounting for low quality reads and index hopping
        *input: 4 multiplexed FASTQ files
        *output: 1 file for each read + index 1 + index 2"""

    if index 1 in R2 contains a converted_phred mean lower than 20:
            write read 1, index 1, index 2 to low_quality_read_1 bin with attached header

    elif index 2 in R3 contains a converted_phred mean lower than 20:
            write read 2, index 1, index 2 to low_quality_read_2 bin with attached header 

    elif index sequence 1 in R2 contains "N":
            write read 1, index 1, index 2 to undetermined_read_1 bin with attached header

    elif index sequence 2 in R3 contains "N":
            write read 1, index 1, index 2 to undetermined_read_2 bin with attached header

    elif index 1 == index 2:
            write read 1, index 1, index 2 to read1_good bin with attached header
            write read 2, index 1, index 2 to read2_good bin with attached header
    return files listed above 
