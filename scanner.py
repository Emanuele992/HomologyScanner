'''
The get_bases_wIndex function returns a string representing the genomic sequence from the FASTA file around the queried position, where the length of the sequence is specified by n_bases. The sequence is centered around the queried position, with n_bases bases on either side.
'''

def get_bases_wIndex(file_path, chromosome, position, n_bases, fai):
    query = ""
    build = open(file_path, "r")
    fasta_index = open(fai, "r")
    start = 0

    try:
        position = int(position)
    except ValueError:
        raise ValueError("Chromosome position must be an integer")
    
    if "chr" not in chromosome:
        chromosome = "chr" + chromosome
    for line in fasta_index:
        if chromosome == line.strip().split("\t")[0]:
            start = int(line.strip().split("\t")[1])
    position_to_reach = start + (int(position) - n_bases) + int((int(position) - n_bases)/50)
    build.seek(position_to_reach, 0) # go to the position to reach from the beginnig of the file
    query = build.read(n_bases * 2 + int((n_bases * 2)/50))
    build.close()
    fasta_index.close()
    query = query.replace("\n", "")
    query = '\n'.join(query[i:i+50] for i in range(0, len(query), 50))
    return(query)

def get_bases(file_path, chromosome, position, n_bases): 
    query = ""
    chromosome_position = 0

    try:
        position = int(position)
    except ValueError:
        raise ValueError("Chromosome position must be an integer")
    
    if "chr" not in chromosome:
        chromosome = "chr" + chromosome
    with open(file_path, "r") as build:
        line = build.readline()
        if chromosome != "chr1":     
            while line:
                line = build.readline()
                if line.startswith('>'):
                    if line.strip()[1:] == chromosome:
                        break
        chromosome_position = build.tell()    
        position_to_reach = chromosome_position + (int(position) - n_bases) + int((int(position) - n_bases)/50)
        build.seek(position_to_reach, 0)
        query = build.read(n_bases * 2 + int((n_bases * 2)/50))
        query = query.replace("\n", "")
        query = '\n'.join(query[i:i+50] for i in range(0, len(query), 50))
        return(query)

'''
Simple blat launcher with os.system. This Function requires conda env blat activation. 
For blat version see conda env info.
'''

def blat_launcher(reference_fa, input_fa, outputname):

    from os import system

    print("Launching blat!")
    # blat launch
    command_to_run = "blat -noHead " + reference_fa + " " + input_fa + " " + outputname + ".psl -out=psl -tileSize=18"
    system(command_to_run)
    

'''
Function description
'''

def homology_score_calculator(match_list, mismatch_list, window):
    match_mismatch_ratio_list = []
    sum_ratios = 0.0
    for i in range(0,len(match_list)):
        if int(match_list[i]) <= 0.90 * 2 * window:
            continue
        if int(mismatch_list[i]) == 0:
           mismatch_list[i] = 1
        match_mismatch_ratio_list.append(float(int(match_list[i])/int(mismatch_list[i])))
    if len(match_mismatch_ratio_list) == 0:
        return(0)
    for ratio in match_mismatch_ratio_list:
         sum_ratios = sum_ratios + ratio
    score = (int(window) * 2)/(sum_ratios/len(match_list))
    return(score)

'''
Function description
'''

def get_blat_info(blat_out):

    match_index = 0
    mismatch_index = 1
    position_index = 9
    chromosomes_index = 13
    position_dict = {}

    with open(blat_out + ".pls", "r") as blat_out:
        for position in blat_out:
            position = position.strip().split("\t")
            if position[position_index] not in position_dict.keys():
                position_dict[position[position_index]] = []
                position_dict[position[position_index]].append((position[match_index], position[mismatch_index], position[chromosomes_index]))
            else:
                position_dict[position[position_index]].append((position[match_index], position[mismatch_index], position[chromosomes_index]))

    return(position_dict)


'''
fastq generator from:
TODO PASTE GIT URL HERE TODO
'''

def fastq_gen(fasta, n_reads, window, temporary_folder):
    from os import system, path
    print("FASTA pre-processing...")
    
    if path.isfile(fasta):
        input_fa = open(fasta, "r")
    else:
        input_fa = fasta

    with open(temporary_folder + "/processed.fa", "w") as out_fa:
        for line in input_fa:
            if line.startswith(">"):
                out_fa.write(line)
                continue
            out_fa.write(line.strip())
    
    if path.isfile(fasta):
        input_fa.close()

    print("Generating fastq...")
    command = "python fastq_generator.py generate_mapped_fastq_SE " + temporary_folder + "/processed.fa " + str(window) + " " + str(n_reads) + " > " + temporary_folder + "/fake.fastq"
    system(command)
    print("Generated!")
    print("Mutating fastq...")
    command = "python mutation_generator.py " + temporary_folder + "/fake.fastq " + temporary_folder + "/fake_mutated.fastq " + str(window)
    system(command)
    command = "rm " + temporary_folder + "/fake.fastq"
    system(command) 

'''
MAPPER
'''

def mapper(bwa_path: str, t: int, temporary_folder: str) -> None:
    import os

    fastq = temporary_folder + "/fake_mutated.fastq.gz"
    if not os.path.isfile(fastq):
        print("gzipping...")
        os.system(f"gzip {fastq[:-3]}")
    
    bam_file = f"{fastq[:-3]}.bam"
    sorted_bam = f"{fastq[:-3]}_sorted.bam"
    sorted_bai = f"{fastq[:-3]}_sorted.bai"
    sorted_bed = f"{fastq[:-3]}_sorted.bed"

    # bwa mem -M -t 8 $BWA ${fastq}.gz > ${fastq%.*}.bam --> Reference command!
    bwa_command = f"bwa mem -M -t {t} {bwa_path} {fastq} > {bam_file}"
    print("Launching bwa-mem...")
    os.system(bwa_command)
    print("Done!")

    # picard SortSam I=${fastq%.*}.bam O=${fastq%.*}_sorted.bam SORT_ORDER=coordinate --> Reference command!
    sort_command = f"picard SortSam I={bam_file} O={sorted_bam} SORT_ORDER=coordinate"
    print("Sorting...")
    os.system(sort_command)
    print("Done!")

    # picard BuildBamIndex I=${fastq%.*}_sorted.bam O=${fastq%.*}_sorted.bai VALIDATION_STRINGENCY=LENIENT
    index_command = f"picard BuildBamIndex I={sorted_bam} O={sorted_bai} VALIDATION_STRINGENCY=LENIENT"
    print("Creating bam index (for IGV/UCSC visualization)...")
    os.system(index_command)
    print("Done!")

    # bamToBed -i ${fastq%.*}_sorted.bam > ${fastq%.*}_sorted.bed
    print("converting bam to bed...")
    bed_command = f"bamToBed -i {sorted_bam} > {sorted_bed}"
    os.system(bed_command)
    print("Done!")

def read_counter(bed_file, chromosome: str, position: int, window: int) -> None:

    def mean(numbers):
            if (len(numbers) == 0):
                return "NA"
            total = 0
            for x in numbers:
                total += x

            return round(total/len(numbers), 2)


    with open(bed_file, "r") as bed:
        results = []
        count_mapped = 0
        count_mismapped = 0
        mapped_qual = []
        mismapped_qual = []
        for line in bed:
            coordinates = line.strip().split("\t")[0:3]
            mapping_quality = line.strip().split("\t")[-2]
            if coordinates[0] == chromosome and position in range(int(coordinates[1]), int(coordinates[2])):
                count_mapped += 1
                mapped_qual.append(int(mapping_quality))
            else:
                count_mismapped += 1
                mismapped_qual.append(int(mapping_quality))
        results = [count_mapped, count_mismapped, mean(mapped_qual), mean(mismapped_qual), round((count_mismapped/(count_mapped + count_mismapped))*100, 2)]
        return(results)
    

def modify_base(ref_fasta, query, position, base, out_file):
    chromosome_position = 0
    out = open(out_file, "w")
    with open(ref_fasta, "rb") as build:
        line = build.readline().decode()
        if line != query: # avoid skipping first line     
            while line:
                line = build.readline().decode()
                if line == query:
                    break
        chromosome_position = build.tell()
        position_to_reach = chromosome_position + int(position) + int(int(position)/50) # correct the position counting "\n"s 1 every 50 bases
        build.seek(0,0)
        string = build.read(position_to_reach - len(base))
        if build.read(1).decode() == "\n":
            build.seek(-2,1)
            out.write(string.decode()[:-len(base)] + base)
            build.read(1)
        else:
            out.write(string.decode() + base) 
        for line in build.readlines():
            out.write(line.decode())
    out.close()


def validate_sequence(sequence):
    from re import match
    if not match('^[ATGC]+$', sequence) or len(sequence) > 1:
        raise ValueError("The provided alternative base is invalid. Only single-character bases A, T, C, and G are allowed as alternatives.")