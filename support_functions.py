"""
Simple spinner object (from stackoverflow)
"""
import sys
import time
import threading

class Spinner:   
    
    busy = False
    delay = 0.1

    @staticmethod
    def spinning_cursor():
        while 1: 
            for cursor in '|/-\\': yield cursor

    def __init__(self, delay=None):
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay): self.delay = delay

    def spinner_task(self):
        while self.busy:
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()

    def __enter__(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()

    def __exit__(self, exception, value, tb):
        self.busy = False
        time.sleep(self.delay)
        if exception is not None:
            return False


def move_to_endline(file):
    while True:
        if file.read(1) == "\n":
            return(file.tell())
        else:
            continue

"""
fa chromosome indexer 
"""

def fa_chr_index(file):
    fai = open(file + ".fai", "w")
    with open(file, "r") as f:
        line = f.readline()
        while line:
            if line.startswith(">"):
                fai.write(line.strip().split(">")[1] + "\t" + str(f.tell()) + "\n")
            line = f.readline()

"""
make ooc: deleterious with short mapping


def make_ooc(file):
    
    from os import system
    system("blat -makeOoc=11.ooc " + file) 

"""

"""
    Given a genomic position in the format "chr1 1234" and a path to a GTF file, returns the distance to the nearest
    exon or the distances of the position to the exon limits if it falls within an exon.

    Args:
        genomic_position_str (str): The genomic position to search for in the GTF file in the format "chr1 1234".
        gtf_path (str): The path to the GTF file.

    Returns:
        int: The distance to the nearest exon or the distances of the position to the exon limits if it falls
        within an exon.
"""

def get_distance_to_nearest_exon(genomic_position_str, gtf_path):
    
    # Parse the chromosome and position from the genomic position string
    chromosome, position = genomic_position_str.split()
    position = int(position)
    min_distance_absolute = [10000000000,'+']
    # Open the GTF file
    with open(gtf_path) as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            # Parse the GTF line
            parts = line.strip().split('\t')
            feature_type = parts[2]
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
           
            

            # Check if the feature type is an exon and if the chromosome matches
            if feature_type == 'exon' and chrom == chromosome:
                # Check if the position is within the exon
                if start <= position <= end:
                    # Return the distance to the exon limit
                    return [abs(position - start), abs(position - end), strand]

                # If the position is not within the exon, calculate the distance to the nearest exon
                distance_to_start = abs(position - start)
                distance_to_end = abs(position - end)
                min_distance = min(distance_to_start, distance_to_end)
                
                if min_distance < min_distance_absolute[0]:
                    min_distance_absolute[0] = min_distance
                    min_distance_absolute[1] = strand
                
                # Return the distance to the nearest exon
            #if not line:
    return min_distance_absolute