# coding=utf-8

from random import randint
import sys

def generate_mutated_fastq(file_path: str, output_path: str, n_bases: int) -> None:
    mutations = []
    N_mutated_sites = n_bases/33 # respect the 3 mutation every 100 bp
    N_mutated_sites = int(N_mutated_sites) # remove . digits
    bases = ["A", "T", "C", "G"]
    with open(output_path, "w") as fastq_mutated:
        with open(file_path, "r+") as fastq:
            for line in fastq:
                if not line.startswith("@") and not line.startswith("+") and not line.startswith("I"):
                    mutations = []
                    for _ in range(N_mutated_sites):
                        mutations.append(randint(0, n_bases - 1)) # this cycle creates up to "N" mutations sites
                    for mutation in mutations:
                        base_pos = randint(0, 3)
                        base = bases[base_pos]
                        line = line[:mutation] + base + line[mutation+1:]
                    fastq_mutated.write(line)
                else:
                    fastq_mutated.write(line) 

def main():
    file = sys.argv[1]
    file_out = sys.argv[2]
    n_bases = sys.argv[3]
    generate_mutated_fastq(file, file_out, int(n_bases))

if __name__ == "__main__":
    sys.exit(main())