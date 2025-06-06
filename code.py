#importing necessary modules/packages

import sys #for command line arguments 
import re #for obtaining cigar string patterns

'''creating command line to obtain user input of files using sys and
assigning variables to user input files'''

sam_file = sys.argv[1] #1st command line argument taking the sam file
input_table = sys.argv[2] #2nd command line argument taking gene input table

#creating function if the cigar string contains N indicating a split read

def split_reads(cigar):
    for match in re.finditer(r'(\d+)([MIDNS])', cigar): #using regex to find pattern in cigar string
        if match.group(2) == 'N': #checking for skipped regions
            return True #presence of N means the read is split as no matches were found
    return False #no skipped regions means the read is not split

#checking the function with assert statements

assert split_reads('16M10I3M') == False #as no N, read is not split
assert split_reads('20M3I4N') == True #as there is N, read is split

#creating function to get junction positions from cigar that has split reads

def junction(pos, cigar):
    junctions = [] #list of start and end positions of junction
    prev_end = pos  #tracking the previous end position (to add the next bases length to). here, starting at the initial position of gene
    
    for match in re.finditer(r'(\d+)([MIDNS])', cigar): #using regex to find pattern in cigar string
        value = int(match.group(1))  #length of the bases
        letter = match.group(2) #result of aligning bases to the reference
        
        if letter in ['M', 'D']:  #incase of match or deletion (as insertions and soft clipped bases not included)
            prev_end = prev_end + value  #adding length of the match or deletion to previous end position 
        
        elif letter == 'N':  #incase of skipped region i.e. junction
            start_pos = prev_end #junction position starts at where the match/deletion ended 
            end_pos = start_pos + value #junction position ends after the length of skipped region
            junctions.append((start_pos, end_pos)) #appending junction start and end positions to a list as tuples for future use
            prev_end = end_pos  #updating prev_end for next junction if any
    return junctions

#checking the function with assert statements
assert junction(10, '20M3I4N') == [(30, 34)] #for 1 junction in cigar
assert junction(1,'16M10I3M20N5M2N') == [(20,40),(45,47)] #for 2 junctions in cigar
assert junction(1, '1M5I10N2M6D4N5I4M1N') == [(2,12),(20,24),(28,29)] #for 3 junctions in cigar

#creating function to count the reads supporting same junction
def count_reads(junctions): #using list of tuples of junction position
    junctions.sort() #arranging so that reads for same junction are grouped together
    readcount = [] #list for storing the junction and the number of reads supporting it
    count = 1 #starting the count at 1 as there is atleast 1 read initially
    
    for i in range(1, len(junctions)): # obtaining indices for iterating over each junction
        if junctions[i] == junctions[i-1]: #checking if junction is same as the previous one
            count += 1 #if same, adding 1 to the read count
        else:
            readcount.append((junctions[i-1], count)) #if  junction is not same as the previous one, storing it and its readcount in the list
            count = 1 #reseting the count to 1 for the next junction
    
    readcount.append((junctions[-1], count)) #to ensure last read is counted
    return readcount

#checking the function with assert statement
junctions1 = [(20, 16), (47, 94), (20, 16), (32, 40), (47, 94), (20, 16)]
assert count_reads(junctions1) == [((20, 16), 3), ((32, 40), 1), ((47, 94), 2)]

#creating function to separate the 4 elements in location column of gene input table
def location_parse(location): 
    chrom, start_end_strand = location.split(':') #extracting chromosome 
    start, end_strand = start_end_strand.split('..') #extracting start position
    end, strand = end_strand.split('(') #extracting end position and strand
    strand = strand.replace(')','') #removing ) from the strand
    start = int(start.replace(',', '')) #removing any , present between the digits of position
    end = int(end.replace(',', ''))  #removing any , present between the digits of position
    return chrom, start, end, strand

#checking the function with assert statements
assert location_parse('TGME49_chrVIII:10000..15000(+)') == ('TGME49_chrVIII', 10000, 15000, '+') #positions without ,
assert location_parse('TGME49_chrX:1,234,567..1,234,789(-)') == ('TGME49_chrX', 1234567, 1234789, '-') #positions with ,

#processing the SAM file
try: #opening if the file exists
    with open(sam_file) as sam:
        junctions = [] #to store the junction positions as tuples
        junc_chrom = [] #to store the chromosome the junctions are on

        for line in sam:
            if line.startswith('@'): #removing headers
                continue

            try: #separating if the columns are in the correct format
                columns = line.strip().split('\t') #separating the columns into a variable
                chrom, pos, cigar, nh = columns[2], int(columns[3]), columns[5], columns[-1]  #assigning variables to the specific columns needed
                
                if nh == 'NH:i:1' and split_reads(cigar): #obtaining split reads aligned only once using split_reads function
                    junctions.append(junction(pos, cigar)) #obtaining the junction positions using junctions function and storing it as tuples in a list
                    junc_chrom.append(chrom) #storing the chromosome associated with the junction in a list

            except ValueError:
                print(f'Error processing line in {sam_file}')
                continue #skipping the line in case of error in formatting and informing the user

        readcount = count_reads(junctions)  #counting the number of reads supporting the function
        
except FileNotFoundError:
    print(f'Error: The file {sam_file} was not found.')
    sys.exit(1) #exiting if the file doesn't exist while informing the user

#processing the gene input file
try: #opening if the file exists
    with open(input_table) as input:
        with open('3027696P.txt', 'w') as output: #creating an output file with student id
            output.write(f'Gene ID\tJunction Start\tJunction End\tNo. of reads supporting the junction\n') #adding a header line so content of the output file are clear
            chrom_list, start_list, end_list, strand_list = [], [], [], [] #creating empty lists to store the gene information
            prev_gene_id = None  #tracking previous gene id so that when gene changes a new line can be inserted
            
            header = input.readline()  #skipping header
            for line in input: 
                try: #separating if the columns are in the correct format
                    gene_id, transcript_id, location = line.strip().split('\t') #assigning variables to the columns in input table
                    chrom, start, end, strand = location_parse(location) #separating the elements in location column and storing them in respective lists
                    chrom_list.append(chrom)
                    start_list.append(start)
                    end_list.append(end)
                    strand_list.append(strand)

                    if gene_id != prev_gene_id: #checking if the gene id has changed from previous one during iteration
                        output.write("\n") #adding a new line when gene id changes i.e. all its junctions have been written
                        #adding a new line when a gene has 0 junctions as well
                        prev_gene_id = gene_id #resetting the previous gene id to current one to start writing the junctions for the new gene

                    for junc_list, count in readcount: #separating junction positions list and their read count in the tuple into variables
                        for junc in junc_list: #iterating over each tuple in the list
                            junc_start, junc_end = junc  # obtaining start and end positions of the junction and assigning them to variables
                        
                            for i in range(0, len(chrom_list)): #obtaining indices to iterate over each junction position
                                if chrom_list[i] in junc_chrom: #if chromosome in the input table is present in the junction chromosome list, then checking whether that junction is within the boundaries of the gene
                                    #checking whether the start and end position of the junction falls in between the start and end positions of the gene
                                    if (start_list[i] <= int(junc_start) <= end_list[i] and
                                        start_list[i] <= int(junc_end) <= end_list[i]): 
                                        output.write(f"{gene_id}\t{junc_start}\t{junc_end}\t{count}\n") #writing the id, positions and read count of the junctions that are within the boundary to the output file
                
                except ValueError:
                    print(f'Error processing line in {input_table}')
                    continue #skipping the line in case of error in formatting and informing the user

except FileNotFoundError:
    print(f'Error: The file {input_table} was not found.')
    sys.exit(1)  #exiting if the file doesn't exist while informing the user
