'''
Generate's CT file for Circle Plot based on protein FASTA sequence

usage: python generateBaseFile.py apoe.fasta
'''

import sys

input_name = sys.argv[-1]
infile = open(input_name, "r")
protein_name = input_name.split('.')[0]
output_name = protein_name + ".ct"
outfile = open(output_name, "w")

#sequence = "KVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVCGRLVQYRGEVQAMLGQSTEELRVRLASHLRKLRKRLLRDADDLQKCLAVYQAGAREGAERGLSAIRERLGPLVEQGRVRAATVGSLAGQPLQERAQAWGERLRARMEEMGSRTRDRLDEVKEQVAEVRAKLEEQAQQIRLQAEAFQARLKSWFEPLVEDMQRQWAGLVEKVQAAVGTSAAPVPSDNH"
sequence = ""

infile.readline()
for line in infile:
    sequence += line.strip()
infile.close()

length = len(sequence)
header = str(length) + " " + protein_name + " file #\n"
outfile.write(header)

for i in range(length):
    text = str(i+1) + " " + sequence[i] + " " + str(i-1+1) + " " + str(i+1+1) + " 0 " + str(i+1) + "\n"
    outfile.write(text)
