'''
Funzione per ottenere i/il fasta della/e sequenza/e da analizzare con blat
la funzione così com'è non va bene è troppo lenta (circa 22/23 sec per query da 150bp per 37 positioni 13.45 min troppissimo)
vanno provate due soluzioni:
    1. funzione seek di python puntando ai bit pos cromosomica start (query - window) --> pos cromosomica end (query + window)
        f.seek(n bit, from what (0 beginning file, 1 current position, 2 end of file)) 
            --> qualcosa tipo f.seek((position-window) - (position-window)/50 (per non contare i \n), 1 dove 1 è l'inizio della riga dopo aver trovato il chr giusto)
        f.tell() tell where you are (absolutely) not so useful honestly
    2. awk ma non ho ancora idea di come.
'''

def get_bases_wIndex(file_path, chromosome, position, n_bases, fai):
    query = ""
    build = open(file_path, "r")
    fasta_index = open(fai, "r")
    start = 0
    if "chr" not in chromosome:
        chromosome = "chr" + chromosome
    for line in fasta_index:
        if chromosome == line.strip().split("\t")[0]:
            start = int(line.strip().split("\t")[1])
    position_to_reach = start + (int(position) - n_bases) + round((int(position) - n_bases)/50)
    build.seek(position_to_reach, 0) # go to the position to reach from the beginnig of the file
    query = build.read(n_bases * 2 + round((n_bases * 2)/50))
    build.close()
    fasta_index.close()
    return(query)

def get_bases(file_path, chromosome, position, n_bases): 
    query = ""
    chromosome_position = 0
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
        position_to_reach = chromosome_position + (int(position) - n_bases) + round((int(position) - n_bases)/50)
        build.seek(position_to_reach, 0)
        query = build.read(n_bases * 2 + round((n_bases * 2)/50))
        
        return(query)

'''
Simple blat launcher with os.system. This Function requires conda env blat activation. 
For blat version see conda env info.
'''

def blat_launcher(reference_fa, input_fa, outputname):

    from os import system

    print("Launching blat!")
    # blat launch
    command_to_run = "blat -noHead " + reference_fa + " " + input_fa + " " + outputname + ".pls"
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

def fastq_gen(fasta, n_reads, window):
    print("FASTA pre-processing...")
    
    with open("processed.fasta", "w") as out_fa:
        with open(fasta, "r") as input_fa:
            for line in input_fa:
                if line.startswith(">"):
                    out_fa.write(line)
                    continue
                out_fa.write(line.strip())

    print("Generating fastq...")
    command = "python fastq_generator.py generate_mapped_fastq_SE processed.fasta " + str(window) + " " + str(n_reads) + " > fake.fastq"
    print("Generated!")
    print("Mutating fastq...")
    command = "python mutation_generator.py fake.fastq fake_mutated.fastq " + str(window)
    from os import system
    system(command)
    #command = "rm fake.fastq"


