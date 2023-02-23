 #!/Usr/bin/env python

import sys, argparse, scanner, support_functions
from time import time
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 


parser = argparse.ArgumentParser(
    description="TODO"
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
    help="tab-separated file with the position to check, one position per line example: \"1 1234\", if provided check the homology score for all the position in the file."
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
    type=str,
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

'''
# DEV mode

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

'   /$$   /$$                                   /$$                                      /$$$$$$                                                             
'  | $$  | $$                                  | $$                                     /$$__  $$                                                            
'  | $$  | $$  /$$$$$$  /$$$$$$/$$$$   /$$$$$$ | $$  /$$$$$$   /$$$$$$  /$$   /$$      | $$  \__/  /$$$$$$$  /$$$$$$  /$$$$$$$  /$$$$$$$   /$$$$$$   /$$$$$$ 
'  | $$$$$$$$ /$$__  $$| $$_  $$_  $$ /$$__  $$| $$ /$$__  $$ /$$__  $$| $$  | $$      |  $$$$$$  /$$_____/ |____  $$| $$__  $$| $$__  $$ /$$__  $$ /$$__  $$
'  | $$__  $$| $$  \ $$| $$ \ $$ \ $$| $$  \ $$| $$| $$  \ $$| $$  \ $$| $$  | $$       \____  $$| $$        /$$$$$$$| $$  \ $$| $$  \ $$| $$$$$$$$| $$  \__/
'  | $$  | $$| $$  | $$| $$ | $$ | $$| $$  | $$| $$| $$  | $$| $$  | $$| $$  | $$       /$$  \ $$| $$       /$$__  $$| $$  | $$| $$  | $$| $$_____/| $$      
'  | $$  | $$|  $$$$$$/| $$ | $$ | $$|  $$$$$$/| $$|  $$$$$$/|  $$$$$$$|  $$$$$$$      |  $$$$$$/|  $$$$$$$|  $$$$$$$| $$  | $$| $$  | $$|  $$$$$$$| $$      
'  |__/  |__/ \______/ |__/ |__/ |__/ \______/ |__/ \______/  \____  $$ \____  $$       \______/  \_______/ \_______/|__/  |__/|__/  |__/ \_______/|__/      
'                                                             /$$  \ $$ /$$  | $$                                                                            
'                                                            |  $$$$$$/|  $$$$$$/                                                                            
'                                                             \______/  \______/  

================================================================================================================================================"""

def main():

    print(welcome_message)

    if args.fasta_index_generator == True:
        print("Fasta index generator mode, this mode disable all other parameters!")
        with support_functions.Spinner():
            support_functions.fa_chr_index(args.reference)
        print ("Index generated! Exiting...")
        sys.exit()
    
    '''
    # DEV mode

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
    
    
    input_fa = "input.fa"

    # -f and -c/-p
    if args.input_file is None and args.chromosome is None and args.position is None:
        parser.error("the following arguments are required: -f/--input_file or -c/--chromosome and -p/--position")

    # -c and -p
    if args.input_file is None and args.chromosome is not None and args.position is not None:
        with support_functions.Spinner():
            query = scanner.get_bases(args.reference, args.chromosome, args.position, args.window)
        
        with open(input_fa, "w") as file_to_query:
            file_to_query.write(">query_" + args.chromosome + "_" + args.position + "_" + str(args.window) + "\n" + query)
        
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
                    else:
                        match_list.append(homology_region[0])
                        mismatch_list.append(homology_region[1])
                    
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    score = scanner.homology_score_calculator(match_list, mismatch_list, args.window)
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + str(score) + "\n")
                    match_list = []
                    mismatch_list = []
                else:
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t0\n")
                
                print("Done!")


    # -f
    if args.input_file is not None and args.chromosome is None and args.position is None:
        file_to_query = open(input_fa, "w")
        print("file provided! message to put here!")
        
        with open(args.input_file) as position_list:
            for line in position_list:
                split_line = line.strip().split("\t")
                chromosome =  split_line[0]
                position = split_line[1]
                print ("query: " + chromosome + " " + position, end = "... ")
                file_to_query.write(">query_" + chromosome + "_" + position + "_" + str(args.window) + "\n")
                with support_functions.Spinner():
                    if args.fasta_index is None:
                        start = time()
                        query = scanner.get_bases(args.reference, chromosome, position, args.window)
                        end = time()
                    else:
                        start = time()
                        query = scanner.get_bases_wIndex(args.reference, chromosome, position, args.window, args.fasta_index)
                        end = time()
                print("elapsed time: " + str(round(end - start, 2)) + "s")
                file_to_query.write(query + "\n")
                
        file_to_query.close()
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
                    else:
                        match_list.append(homology_region[0])
                        mismatch_list.append(homology_region[1])
                    
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    score = scanner.homology_score_calculator(match_list, mismatch_list, args.window)
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t" + str(score) + "\n")
                    match_list = []
                    mismatch_list = []
                else:
                    print("position (chr,position,window): " + ",".join(position.split("_")[1:len(position.split("_"))]), end = "... ")
                    output_file.write("\t".join(position.split("_")[1:len(position.split("_"))]) + "\t0\n")
                
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

