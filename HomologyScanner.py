#!/Usr/bin/env python

import sys, argparse, scanner, support_functions, os
from time import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 


parser = argparse.ArgumentParser(
    description="This program calculates the homology level of genomic positions on an input sequence (in FASTA format) and a window around the position (which can be assumed to be the length of the reads),\
        generating a score that indicates the probability of read mismapping. The homology level is calculated using the output of BLAT (Conda version). \
        This program also has a simulation mode (currently only available for a single position) that allows you to generate a number 'n' of reads at the requested position (twice the window size) with the same length as the window.\
        Afterwards, these reads are randomly mutated by introducing a number of mutations that respects the following ratio: 100bp/3mut. These newly generated reads are then mapped using BWA-MEM2, and the ratios between the number of correctly mapped reads versus those mismapped or unmapped are calculated."
)

parser.add_argument(
    "-r",
    "--reference",
    type=str,
    default=None,
    help="the reference genome fasta",
    required=True
)

parser.add_argument(
    "-w",
    "--window",
    type=int,
    default=150,
    help="window around the queried base to check homology, default 150bp"
)

parser.add_argument(
    "-f",
    "--input_file",
    type=str,
    default=None,
    help="tab-separated file with the position to check, one position per line example: \"1 1234\", \
        if provided check the homology score for all the position in the file."
)

parser.add_argument(
    "-c",
    "--chromosome",
    type=str,
    default=None,
    help="chromosome to check"
)

parser.add_argument(
    "-p",
    "--position",
    type=int,
    default=None,
    help="position to check"
)

parser.add_argument(
    "-b",
    "--blat_output_prefix",
    type=str,
    default="blat_output",
    help="blat output prefix"
)

parser.add_argument(
    "-o",
    "--output_prefix",
    type=str,
    default="homology_score_output",
    help="homology score output prefix"
)

parser.add_argument(
    "-i",
    "--fasta_index",
    type=str,
    default=None,
    help="fasta index created with fasta index creator function"
)

parser.add_argument(
    "-I",
    "--fasta_index_generator",
    action='store_true',
    default=None,
    help="fasta index chromosome generator"
)

parser.add_argument(
    "-R",
    "--read_generation_mode",
    action='store_true',
    default=None,
    help="If the read generation mode is set to true, the program will generate reads of the specified window length from the queried position (default 1000) and introduce random mutations (3-100 bp) into them. \
        The resulting reads will be mapped (bwa-mem2) onto the reference genome specified by the -r argument, and the program will output the percentage of reads that were not correctly mapped.\
        Enabled only if -c and -p arguments are provided"
)

parser.add_argument(
    "-n_reads",
    "--reads_number",
    type=int,
    default=1000,
    help="Number of read that have to be generated (default 1000)."
)

parser.add_argument(
    "-BWA",
    "--bwamem_index_generator_mode",
    action="store_true",
    default=None,
    help="bwa-mem index creator function, if true the program creates a folder with the index necessary to run bwa-mem required with the 'read generation mode' on."
)

parser.add_argument(
    "-bwa",
    "--bwamem_index",
    type=str,
    default=None,
    help="bwa-mem index path. /path/to/bwa-mem_index/build_folder/build"
)

parser.add_argument(
    "-t",
    "--treads",
    type=int,
    default=4,
    help="number of tread/cpus available. This parameter affects only bwa mem."
)

parser.add_argument(
    "-m",
    "--mutation_position",
    type=int,
    default=None,
    help="Relative position that have to be mutated i.e. queried position: chromosome 1 position 1000, -m 10 --> modified base chromosome 1 position 1010"
)

parser.add_argument(
    "-a",
    "--alternative",
    type=str,
    default=None,
    help="Base or bases that have to be inserted in the specified position, it requires '-m'"
)
'''
# NOT WORKING

parser.add_argument(
    "-ooc",
    "--ooc",
    action='store_true',
    default=None,
    help="Over-occurring 11-mers file"
)

parser.add_argument(
    "-OOC",
    "--ooc_blat_generator",
    action='store_true',
    default=None,
    help="Over-occurring 11-mers file generator"
)

'''

args = parser.parse_args()

welcome_message = """
================================================================================================================================================

Welcome to Homology Scanner! The first program that computes a variant-calling specific homology score!

================================================================================================================================================
"""


endMessage = """
================================================================================================================================================
Thanks for using...

   /$$   /$$                                   /$$                                      /$$$$$$                                                             
  | $$  | $$                                  | $$                                     /$$__  $$                                                            
  | $$  | $$  /$$$$$$  /$$$$$$/$$$$   /$$$$$$ | $$  /$$$$$$   /$$$$$$  /$$   /$$      | $$  \__/  /$$$$$$$  /$$$$$$  /$$$$$$$  /$$$$$$$   /$$$$$$   /$$$$$$ 
  | $$$$$$$$ /$$__  $$| $$_  $$_  $$ /$$__  $$| $$ /$$__  $$ /$$__  $$| $$  | $$      |  $$$$$$  /$$_____/ |____  $$| $$__  $$| $$__  $$ /$$__  $$ /$$__  $$
  | $$__  $$| $$  \ $$| $$ \ $$ \ $$| $$  \ $$| $$| $$  \ $$| $$  \ $$| $$  | $$       \____  $$| $$        /$$$$$$$| $$  \ $$| $$  \ $$| $$$$$$$$| $$  \__/
  | $$  | $$| $$  | $$| $$ | $$ | $$| $$  | $$| $$| $$  | $$| $$  | $$| $$  | $$       /$$  \ $$| $$       /$$__  $$| $$  | $$| $$  | $$| $$_____/| $$      
  | $$  | $$|  $$$$$$/| $$ | $$ | $$|  $$$$$$/| $$|  $$$$$$/|  $$$$$$$|  $$$$$$$      |  $$$$$$/|  $$$$$$$|  $$$$$$$| $$  | $$| $$  | $$|  $$$$$$$| $$      
  |__/  |__/ \______/ |__/ |__/ |__/ \______/ |__/ \______/  \____  $$ \____  $$       \______/  \_______/ \_______/|__/  |__/|__/  |__/ \_______/|__/      
                                                             /$$  \ $$ /$$  | $$                                                                            
                                                            |  $$$$$$/|  $$$$$$/                                                                            
                                                             \______/  \______/  

================================================================================================================================================"""

def main():

    print(welcome_message)

    if args.fasta_index_generator == True:
        print("Fasta index generator mode, this mode disable all other parameters!")
        with support_functions.Spinner():
            support_functions.fa_chr_index(args.reference)
        print ("Index generated! Exiting...")
        sys.exit()

    if args.bwamem_index_generator_mode == True:
        print("BWA-MEM index generator mode, this mode disable all other parameters!")
        print("ATTENTION BWA-MEM INDEX REQUIRES A LARGE AMOUNT RAM (FOR HG19 AT LEAST 18 GB) AND FREE DISK SPACE (12-15 GB)")
        with support_functions.Spinner():
            support_functions.bwa_mem_index(args.reference)
        print ("Index generated! Exiting...")
        sys.exit()
    
    if args.window < 33:
        raise ValueError("Please ensure that the window size is a minimum of 33 bases.") # error to have at least a mutation per read generated

    '''
    # NOT WORKING

    if args.ooc_blat_generator == True:
        print("over-occurring 11-mers generator mode, this mode disable all other parameters!")
        with support_functions.Spinner():
            support_functions.make_ooc(args.reference)
        print ("Ooc file generated! Exiting...")
        from sys import exit
        exit()
    if args.ooc == False:
        parser.error("he following arguments are required: -ooc/--ooc")

    '''
    # support variables
    temporary_folder_path = "./temporary"
    input_fa = temporary_folder_path + "/input.fa"

    if not os.path.exists(temporary_folder_path):
        os.makedirs(temporary_folder_path)
        print("Temporary folder created!")
    else:
        print("Temporary folder already exists... ")

    # -f and -c/-p
    if args.input_file is None and args.chromosome is None and args.position is None:
        parser.error("the following arguments are required: -f/--input_file or -c/--chromosome and -p/--position")

    # -c and -p
    if args.input_file is None and args.chromosome is not None and args.position is not None:
        with open(input_fa, "w") as file_to_query:
            with support_functions.Spinner():
                if args.fasta_index is None:
                    start = time()
                    query = scanner.get_bases(args.reference, args.chromosome, args.position, args.window)
                    end = time()
                else:
                    start = time()
                    query = scanner.get_bases_wIndex(args.reference, args.chromosome, args.position, args.window, args.fasta_index)
                    end = time()
            print("elapsed time: " + str(round(end - start, 2)) + "s")
            file_to_query.write(">query_" + args.chromosome + "_" + str(args.position) + "_" + str(args.window) + "\n" + query)
        
        # mutation generator
        if args.mutation_position is not None and args.alternative is not None:
            modified_input_fa = "mutated_input.fa"
            scanner.validate_sequence(args.alternative)
            query_string = ">query_" + args.chromosome + "_" + str(args.position) + "_" + str(args.window) + "\n" # keep \n it is not stripped in the scanner funciton!
            scanner.modify_base(input_fa, query_string, args.mutation_position, args.alternative, modified_input_fa)
            input_fa = modified_input_fa
        elif (args.mutation_position is not None and args.alternative is None) or (args.mutation_position is None and args.alternative is not None):
            parser.error("please provide both -m/--mutation_position and -a/--alternative parameters.")

        if args.read_generation_mode == True:
            print("Read Generator mode...")
            with support_functions.Spinner():
                scanner.fastq_gen(input_fa, args.reads_number, args.window, temporary_folder_path) 

            print("...fastq generated!")
            if args.bwamem_index is None:
                print("bwa-mem index missing... assuming reference path...")
                index_folder = args.reference.split("/")[-1]
                index_folder = index_folder.split(".")[0]
                bwa =  args.reference + "/" + index_folder
                print(bwa)
                sys.exit()

            with support_functions.Spinner():
                scanner.mapper(args.bwamem_index, args.treads, temporary_folder_path)

            sequencing_results = []
            with support_functions.Spinner():
                sequencing_results = scanner.read_counter(temporary_folder_path + "/fake_mutated.fastq_sorted.bed", args.chromosome, args.position, args.window)

        # Launching BLAT
        scanner.blat_launcher(args.reference, input_fa, args.blat_output_prefix)

        # score
        with support_functions.Spinner():
            position_dict = scanner.get_blat_info(args.blat_output_prefix)
        
        with open(args.output_prefix + ".txt", "w") as output_file:
            match_list = []
            mismatch_list = []

            print("Homology score computation...")

            for position in position_dict.keys():
                if len(position_dict[position]) > 1:
                    for homology_region in position_dict[position]:
                        if homology_region[0] == args.window * 2 and homology_region[1] == "0": # skip queried position
                            continue
                        elif "_" in homology_region[2]: # skip special chromosomes e.g. chrUn_gl000212
                            continue
                    else:
                        match_list.append(homology_region[0])
                        mismatch_list.append(homology_region[1])
                    
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    score = scanner.homology_score_calculator(match_list, mismatch_list, args.window)
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + str(score) + "\t" + "\t".join(map(str, sequencing_results)) + "\n")
                    match_list = []
                    mismatch_list = []
                else:
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + "0" + "\t" + "\t".join(map(str, sequencing_results)) + "\n")
                
                print("Done!")


    # -f
    if args.input_file is not None and args.chromosome is None and args.position is None:
        file_to_query = open(input_fa, "w")
        sequencing_results_dict = {}
        print("file provided!")
        
        mutation_flag = False # this flag check the number of columns of the file, if two flag is False and the mutation generation is skipped

        with open(args.input_file) as position_list:
            for line in position_list:
                split_line = line.strip().split("\t")
                if len(split_line) == 4:
                    mutation_flag = True
                if len(split_line) != 2 and len(split_line) != 4:
                    raise ValueError("Columns number must be 2 or 4")
                chromosome =  split_line[0]
                position = split_line[1]
                query_string = ">query_" + chromosome + "_" + str(position) + "_" + str(args.window) + "\n"
                print("query: " + chromosome + " " + position, end = "... ")
                file_to_query.write(query_string)
                if args.fasta_index is None:
                    start = time()
                    query = scanner.get_bases(args.reference, chromosome, position, args.window)
                    end = time()
                else:
                    start = time()
                    query = scanner.get_bases_wIndex(args.reference, chromosome, position, args.window, args.fasta_index)
                    end = time()
                print("elapsed time: " + str(round(end - start, 2)) + "s")
                if mutation_flag:
                    # Mutating!
                    temporary_query_file_input = temporary_folder_path + "/temporary_query.fa"
                    temporary_query_file_output = temporary_folder_path + "/temporary_query_mutated.fa"
                    mut_pos = split_line[2]
                    try:
                        mut_pos = int(mut_pos)
                    except ValueError:
                        raise ValueError("Mutation position must be an integer")
                    alt = split_line[3] 
                    scanner.validate_sequence(alt) # check if the mutation provided contains inconsistencies
                    
                     # keep \n it is not stripped in the scanner funciton!

                    # temporary output
                    with open(temporary_query_file_input, "w") as temp_query:
                        temp_query.write(query_string)
                        temp_query.write(query + "\n")

                    scanner.modify_base(temporary_query_file_input, query_string, int(mut_pos), alt, temporary_query_file_output)
                    
                    with open(temporary_query_file_output, "r") as temp_query_out:
                        modified_query = ""
                        for line in temp_query_out:
                            if line.startswith(">"): 
                                continue
                            modified_query = modified_query + line
                    query = modified_query

                file_to_query.write(query + "\n")
            file_to_query.close()
            '''
                    if args.read_generation_mode == True:
                        print("Read Generator mode...")

                        # Splitting input fasta in different temporary fasta files

                        with support_functions.Spinner():
                            scanner.fastq_gen(input_fa, args.reads_number, args.window, temporary_folder_path) 

                        print("...fastq generated!")
                        if args.bwamem_index is None:
                            print("bwa-mem index missing... assuming reference path...")
                            index_folder = args.reference.split("/")[-1]
                            index_folder = index_folder.split(".")[0]
                            bwa =  args.reference + "/" + index_folder
                            print(bwa)
                            sys.exit()

                        with support_functions.Spinner():
                            scanner.mapper(args.bwamem_index, args.treads, temporary_folder_path)

                        key = chromosome + "_" + position + "_" + args.window
                        
                        # TO MODIFY
                        if key in sequencing_results_dict.keys():
                            raise KeyError("Ther variant " + " ".join(key.split("_")) + " is present multiple times")

                        with support_functions.Spinner():
                            sequencing_results_dict[key] = scanner.read_counter(temporary_folder_path + "/fake_mutated.fastq_sorted.bed", chromosome, position, args.window)'''
                
        # Launching BLAT
        scanner.blat_launcher(args.reference, input_fa, args.blat_output_prefix) 

        # score
        with support_functions.Spinner():
            position_dict = scanner.get_blat_info(args.blat_output_prefix)

        with open(args.output_prefix + ".txt", "w") as output_file:
            match_list = []
            mismatch_list = []

            print("Homology score computation...")

            for position in position_dict.keys():
                if len(position_dict[position]) > 1:
                    for homology_region in position_dict[position]:
                        if homology_region[0] == args.window * 2 and homology_region[1] == "0": # skip queried position
                            continue
                        elif "_" in homology_region[2]: # skip special chromosomes e.g. chrUn_gl000212
                            continue
                    else:
                        match_list.append(homology_region[0])
                        mismatch_list.append(homology_region[1])
                    
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    score = scanner.homology_score_calculator(match_list, mismatch_list, args.window)
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + str(score)  + "\t" + "\t".join(map(str, sequencing_results_dict[position])) + "\n")
                    match_list = []
                    mismatch_list = []
                else:
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + "0" + "\t" + "\t".join(map(str, sequencing_results_dict[position])) + "\n")
                
                print("Done!")
                

        
    # -c and not -p (and viceversa)
    if (args.chromosome is None and args.position is not None) or (args.chromosome is not None and args.position is None):
        parser.error("the following arguments are required: -f/--input_file or -c/--chromosome and -p/--position")

    # not (-f -c -p)   
    if args.input_file is not None and args.chromosome is not None and args.position is not None:
        parser.error("incopatible arguments, provide only -f/--input_file or -c/--chromosome and -p/--position")

    print(endMessage)

if __name__ == "__main__":
    sys.exit(main())

