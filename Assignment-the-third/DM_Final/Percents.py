#!/usr/bin/env python
import os 

dirname = '/home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-third'
ext = '.fq'

#keys = filename, values = number of records in that file
final_counts_dict = {}

#keys = filename, values = % of the total
percents_dict = {}     
for filename in os.listdir(dirname):      #goes through every file in dirname
    if filename.endswith(ext):            #targets filenames ending with ext
        with open(filename) as f:          
            count = sum(0.25 for line in f)    #easy way to count number of records based on line
            final_counts_dict[filename] = int(count)
#I add a new key:value pair, where I sum the values of both R1 and R2 hopped files
final_counts_dict['Hopped.R1.R2'] = final_counts_dict['Hopped.R1.fq'] + final_counts_dict['Hopped.R2.fq']

#I delete individual hopped files because I was double-counting them when finding the percents
del final_counts_dict['Hopped.R1.fq']
del final_counts_dict['Hopped.R2.fq']

#values = filename, key = % of total reads in each sample
for filename, count in final_counts_dict.items():
    percents_dict[filename] = (count/(sum(final_counts_dict.values()))*100)

#write percents to new file
for filename in os.listdir(dirname):
    if filename.endswith(ext):
        final_report_fh = open(f'Final_Report.txt', 'w')
        final_report_fh.write(f'Filename\tPercent\n')
        for filename, percent in percents_dict.items():
            final_report_fh.write(f'{filename}\t{round(percent,2)}%\n')       #probably don't need 5 decimal points for this