import sys
import sequences
#Use python epitopescan.py seq_file epitope_file output


#generate a dictionary of sequences from a fasta file
infh1 = open(sys.argv[1])
readfasta= sequences.FastaIterator(infh1)

seq_dico={}
seq = readfasta.next()
while seq:
    seq_dico[seq.name]= seq.translation()
    seq=readfasta.next()

infh1.close()

#read a list of epitopes from a file
infh2 = open(sys.argv[2])
epitopes = []
line = infh2.readline()
line = infh2.readline()
while line:
    if len(line) < 2:
	line = infh2.readline()
    line = line.replace('\n', '')
    line = line.replace('\r', '')
    a = line.split(',')
    epitopes.append(a)
    line = infh2.readline()
infh2.close()

patients = seq_dico.keys()


#Checks if each epitope from the list is in a given sequence and generates a binary table of presence/absence of epitopes
out = ''
for p in patients:
    name = p
    out += name
    for e in epitopes:
	    if e[0] in seq_dico[p]:
	        out +='\t1'
	    else:
	        out +='\t0'
    out+='\n'


head = 'Sequence'
for e in epitopes:
    head+= '\t'+e[0]+'_'+e[1]
head +='\n'

out = head + out

#Writes the table into an output file
outfh = open(sys.argv[3],'w')
outfh.write(out)

outfh.close()
