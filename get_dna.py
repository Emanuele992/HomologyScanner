from support_functions import get_sequence_from_genomic_coordinates
import sys

def main():
    chromosome = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    build = sys.argv[4]
    query_fasta = get_sequence_from_genomic_coordinates(chromosome, start, end, build)
    print(query_fasta)

if __name__ == "__main__":
    sys.exit(main())