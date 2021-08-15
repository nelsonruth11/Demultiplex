if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

def validate_base_seq(seq: str,RNAflag: bool =False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    DNAbases = set('ATGCatcg')
    RNAbases = set('AUGCaucg')
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def validate_DNA_seq(DNA): 
    '''See if string is DNA statement'''
    DNA = DNA.lower()
    A = DNA.count('a')
    T = DNA.count('t')
    C = DNA.count('c')
    G = DNA.count('g')
    sum_agtc = A + C + G + T
    
    return sum_agtc == len(DNA)

def gc_content(DNA):
    '''Calculate GC%'''
    DNA.upper()
    GC_tot = DNA.count('G') + DNA.count('C')
    GC_content = GC_tot/(len(DNA))

    return GC_content

def convert_phred(letter):
    """Converts a single character into a phred score"""
    letter = ord(letter) - 33
    return letter

def get_args():
    """Basic args formatting for various occasions. Need args = get_args()\n input_file = args.input_filename...after"""
    parser = argparse.ArgumentParser(description="A program to normalize kmer data")
    parser.add_argument("-f", "--input_filename", help="your filename", required = True)
    parser.add_argument("-t", "--table", help="ENSEMBL table", required = True)
    parser.add_argument("-l", "--length", help="length of read", required = True, type=int)
    parser.add_argument("-k", "--kmer", help="length of kmer", type=int)
    parser.add_argument("-c", "--coverage_limit" ,help="specify coverage", required = True, type=int)
    parser.add_argument("-o", "--output_filename", help="your output filename")

    return parser.parse_args()

def biomart_to_two_dicts(table):
    """Takes biomart table as input. Returns two dicts:
        *prot_id : gene_id 
        *protID : gene_name
        *requires Biomart table to be order: Gene ID, Protein ID, Gene Name
    """
    protID_geneID_dict ={}    #initialize empty dictionary with keys as protein ID and values as geneID
    protID_gene_name_dict = {}     #initialize empty dictionary with keys as protein ID and values as gene name

    with open (table, "r") as fh:     #opens predownloaded biomart table
        while True:
            line = fh.readline().strip()     #strips each line
            if line =='':
                break
            parts = line.split("\t")     #parts is entire header split by tab 
            if len(parts) == 3:
                geneID, protID, gene_name = parts
            if len(parts) == 2:
                geneID, protID = parts
                gene_name = ""                
                assert "ENS" in geneID      #gives assert error if no geneID
                assert "ENS" in protID      #gives assert error if no protID 
            protID_geneID_dict[protID] = geneID
            protID_gene_name_dict[protID] = gene_name

    return protID_geneID_dict, protID_gene_name_dict

def convert_to_two_line_fasta(input_file):
    """Converts multilined FASTA files to two-lined"""
    seq = ""     #initialize seq to ""
    first_line=True
    #gene_id_regex = "ENSDARG[0-9]+"
    head_list = []     #initialize empty list for full headers
    seq_list = []       #initialize list for protein sequences only

    with open(input_file, "r") as fh:     #must be fasta file
        while True:                       
            line = fh.readline().strip()
            if line =='':
                break
            if line[0] == '>':
                if first_line:
                    header=line
                    first_line=False
                else:
                    head_list.append(header)     #appends header 
                    seq_list.append(seq)        #appends seq
                    header = line
                    seq = ""
            else:
                seq += line
    head_list.append(header) #appends the final header line
    seq_list.append(seq)    #appends the final sequence line
    return head_list, seq_list


def rev_comp(seq):
    '''Reverse complements a DNA sequence'''
    seq = seq.upper()
    complement = str()
    for base in seq:
        if base == "A":
            complement = complement + str("T")
        if base == "T":
            complement = complement + str("A")
        if base ==  "C":
            complement = complement + str("G")
        if base == "G":
            complement = complement + str("C")
        if base == "N":
            complement = complement + str("N")
    reverse_complement = complement[::-1]
    return reverse_complement