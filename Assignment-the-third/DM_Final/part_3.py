#!/usr/bin/env python

#import functions
import gzip, argparse
from Bioinfo import rev_comp, convert_phred

#hardcoded source files
# r1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"      #read 1
# r2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"      #index 1
# r3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"      #index 2
# r4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"      #read 2
# indexes = "/projects/bgmp/shared/2017_sequencing/indexes.txt"                        #index file

def get_args():
    parser = argparse.ArgumentParser(description="A program demultiplex paired-end sequencing data")
    parser.add_argument("-f1", "--input_r1", help="read1 filename", required = True)
    parser.add_argument("-f2", "--input_r2", help="index1 filename", required = True)
    parser.add_argument("-f3", "--input_r3", help="index2 filename", required = True)
    parser.add_argument("-f4", "--input_r4", help="read2 filename", required = True)
    parser.add_argument("-i", "--input_index", help="indexes filename", required = True)

    return parser.parse_args()

args = get_args()

input_r1 = args.input_r1
input_r2 = args.input_r2
input_r3 = args.input_r3
input_r4 = args.input_r4
input_i = args.input_index

#create dictionary from indexes.txt file
#key = index name, value = index sequence
index_dict = {}

#create dictionary for output files
#keys = index sequence, values = filehandles for writing
output_fh_dict = {}

############################################
#populate index_dict from .txt file

with open(input_i, "r") as fh:            
    next(fh)       #skips header line
    for line in fh:
        line = line.split()
        index_dict[line[4]] = line[3]

#populate output file dict
for index_sequence, index_name in index_dict.items():
    # eg.: GTAGCGTA.R1.fq for writing
    read1_fh = open(f'{index_sequence}.R1.fq', "w")    
    read2_fh = open(f'{index_sequence}.R2.fq', "w")
    output_fh_dict[index_sequence] = [read1_fh, read2_fh]

#hardcoded these because they'll be used regardless of project
read_1_hopped_fh = open(f'Hopped.R1.fq', 'w')
read_2_hopped_fh = open(f'Hopped.R2.fq', 'w')
read_1_low_quality_fh = open(f'Low_Quality.R1.fq', 'w')
read_2_low_quality_fh = open(f'Low_Quality.R2.fq', 'w')
final_report_fh = open(f'Final_Report.txt', 'w')

################
#creating FASTQ entries with index sequences attached

with gzip.open(input_r1, "r") as fr1, gzip.open(input_r2, "r") as fr2, gzip.open(input_r3, "r") as fr3, gzip.open(input_r4, "r") as fr4:  
    read_1_header_final = ''
    read_1_entry_final = ''

    for read1, index1, index2, read2 in zip(fr1, fr2, fr3, fr4):
        #need to decode lines from compressed files
        read1 = read1.decode('UTF-8')   
        read2 = read2.decode('UTF-8')
        index1 = index1.decode('UTF-8')
        index2 = index2.decode('UTF-8')
        num_lines += 1
        if num_lines % 4 == 1:     #headers only
            read_1_header = read1   #need this step to add index sequences to header
            read_2_header = read2
        if num_lines % 4 == 2:      #sequences only
            read_1_header_final = "-".join([read_1_header.strip(), index1.strip(), rev_comp(index2.strip())])   #adds sequences of indexes to header
            read_1_entry = (f'{read_1_header_final.strip()}\n{read1.strip()}\n+')      #creates FASTQ entry
            read_2_header_final = "-".join([read_2_header.strip(), index1.strip(), rev_comp(index2.strip())])   #adds sequences of indexes to header
            read_2_entry = (f'{read_2_header_final.strip()}\n{read2.strip()}\n+')
        if num_lines % 4 == 0:     #quality lines only
            read_1_entry_final = (f'{read_1_entry}\n{read1.strip()}\n')      #final Read 1 FASTQ entry, to be written later
            read_2_entry_final = (f'{read_2_entry}\n{read2.strip()}\n')      #final Read 2 FASTQ entry, to be written later
            read_1_header_split = read_1_header_final.split("-")
            read_2_header_split = read_2_header_final.split("-")
#############
#filtering by mean Qscore
            #initialize counters for quality scores
            qscore_read_1_sum = 0
            qscore_read_2_sum = 0
            qscore_index_1_sum = 0
            qscore_index_2_sum = 0

            #filtering based on quality of read
            for score in read1.strip():
                qscore_read_1_sum += convert_phred(score)                  
            q_read_1_mean = (qscore_read_1_sum/len(read1))
            if q_read_1_mean < 20:
                read_1_low_quality_fh.write(read_1_entry_final)      #write to R1 low wuality file

                read_1_low_quality_count += 1

            for score in read2.strip():
                qscore_read_2_sum += convert_phred(score)                   
            q_read_2_mean = (qscore_read_2_sum/len(read1))

            if q_read_2_mean < 20:
                read_2_low_quality_fh.write(read_2_entry_final)           #write to R2 low quality file

                read_2_low_quality_count += 1

            #filtering based on quality of index
            for score in index1.strip():
                qscore_index_1_sum += convert_phred(score)          
            q_index_1_mean = (qscore_index_1_sum/len(index1))

            if q_index_1_mean < 30:
                read_1_low_quality_fh.write(read_1_entry_final)              #write to R1 low quality file
                read_2_low_quality_fh.write(read_2_entry_final)              #write to R2 low quality file

            for score in index2.strip():
                qscore_index_2_sum += convert_phred(score)          
            q_index_2_mean = (qscore_index_2_sum/len(index2))

            if q_index_2_mean < 30:
                read_1_low_quality_fh.write(read_1_entry_final)            #write to R1 low quality file
                read_2_low_quality_fh.write(read_2_entry_final)            #write to R2 low quality file

###############
#filtering by if index not in index dictionary or if indexes have hopped
            
            if q_read_1_mean > 20 and q_read_2_mean > 20 and q_index_1_mean > 30 and q_index_2_mean > 30:
                #if indexes are not in the index dictionary, also catches any indexes containing N
                if read_1_header_split[-2] not in index_dict.keys() or read_1_header_split[-1] not in index_dict.keys():
                    read_1_low_quality_fh.write(read_1_entry_final)    #write both entries to low quality file
                    read_2_low_quality_fh.write(read_2_entry_final)
                    read_1_low_quality_count += 1
                    read_2_low_quality_count += 1
                #if indexes are hopped    
                elif read_1_header_split[-2] != read_1_header_split[-1] or read_2_header_split[-2] != read_2_header_split[-1]:   
                    read_1_hopped_fh.write(read_1_entry_final)    #write both entries to hopped file
                    read_2_hopped_fh.write(read_2_entry_final)
                    #print(f'Hopped: {read_1_header_final}'               
                    hopped_count += 1

                #finally we can write to the proper file
                elif read_1_header_split[-2] == read_1_header_split[-1]:
                    output_fh_dict[read_1_header_split[-1]][0].write(read_1_entry_final)     #writes R1 entry to R1 index file
                    read_match_counter += 1
                    output_fh_dict[read_1_header_split[-1]][1].write(read_2_entry_final)     #writes R2 entry to R2 index file
