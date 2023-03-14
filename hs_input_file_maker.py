import sys

chromosome = sys.argv[1]
start_pos = sys.argv[2]
end_pos = sys.argv[3]
filename = sys.argv[4]

with open(filename + ".txt", "w") as output:
    for i in range(int(start_pos),int(end_pos)):
        output.write(chromosome + "\t" + str(i) + "\n")