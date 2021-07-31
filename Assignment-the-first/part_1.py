!/usr/bin/env python

def convert_phred(letter):
    """Converts a single character into a phred score"""
    letter = ord(letter) - 33
    return letter

file = "/home/ndr/bgmp/bioinfo/Bi622/Demultiplex/Assignment-the-first/r1"
#file = "../../test.fastq"
## Use this cell to populate your numpy array of arrays
    
#my_array = np.zeros((101,4000000), dtype=float)
#np.set_printoptions(threshold=np.inf)    #TURN THIS OFF FOR FULL RECORD

num_lines = 0
with open(file, "r") as fh:       #opens and reads file as fh
    for line in fh:               #loops thru file, counter ticks up, line stripped
        num_lines += 1
        line = line.strip()
        if num_lines % 4 == 0: 
            for index, score in enumerate(line):
                my_array[index, (num_lines//4)-1] = Bioinfo.convert_phred(score)
print(my_array)

