import sys
import numpy
from numpy import *
import random
import os
random.seed()

def main():
	motif_length = int(sys.argv[1])
	#print motif_length
	num_variablepositions = int(sys.argv[2])
	#print num_variablepositions
	sequence_length = int(sys.argv[3])
	#print sequence_length
	sequence_count = int(sys.argv[4])
	setnum = int(sys.argv[5])
	path = "dataset"+str(setnum)+"/"
	#print sequence_count
	
	if not os.path.exists(path):
    		os.makedirs(path)
	all_sequences = []
	
	for i in range(sequence_count):
		seq = []
		for k in range(sequence_length):
			j = random.randint(0,3)
			if j==0:
				seq.append("A")
			elif j==1:
				seq.append("C")
			elif j==2:
				seq.append("G")
			elif j==3:
				seq.append("T")
		all_sequences.append(seq)
	for i in range(len(all_sequences)):
		print "".join(all_sequences[i])
	
	motif = []
	motif_map = []
	for i in range(motif_length):
		j = random.randint(0,3)
		if j==0:
			motif.append("A")
			motif_map.append("0")
		elif j==1:
			motif.append("C")
			motif_map.append("0")
		elif j==2:
			motif.append("G")
			motif_map.append("0")
		elif j==3:
			motif.append("T")
			motif_map.append("0")
	
	# create list of num_variablepositions number of *unique* random numbers 
	# in the range [0, motif_length):
	random_positions = random.sample(range(0,motif_length), num_variablepositions) 
	for i in range(num_variablepositions):
		j = random_positions[i]
		motif_map[j]="*";	
	print "".join(motif)
	print "".join(motif_map)
	
	bindingsites = []
	for i in range(sequence_count):
		bs = ""
		for j in range(motif_length):
			if motif_map[j] is "0":
				bs+=motif[j]
			else:	# elif motif_map[j] is "*"
				k = random.randint(0,3)
				if k==0:
					bs+="A"
				elif k==1:
					bs+="C"
				elif k==2:
					bs+="G"
				elif k==3:
					bs+="T"
		bindingsites.append(bs)
	print bindingsites
	locs = []
	for i in range(sequence_count):
		loc = random.randint(0, sequence_length-motif_length)
		for j in range(motif_length):
			all_sequences[i][j+loc]=bindingsites[i][j]
		locs.append(loc)
		
	f = open(path+'sequences.fa','w')
	f2 = open(path+'sites.txt', 'w')
	for i in range(len(all_sequences)):
		print "".join(all_sequences[i])
		label = ">S"+str(i)+"\n";
		sequence = "".join(all_sequences[i])+"\n"
		f.write(label);
		f.write(sequence);
		str2 = "S"+str(i)+"\t"+str(locs[i])+"\n"
		f2.write(str2)
	f.close();
	f2.close()
	
	mainmotif = ""
	for i in range(motif_length):
		if motif_map[i] is "0":
			mainmotif+=motif[i]
		else:
			mainmotif+="*"
	f3 = open(path+'motif.txt', 'w')
	s3text = "MOTIF"+"\t"+str(motif_length)+"\t"+mainmotif
	f3.write(s3text)
	f3.close();
	f4 = open(path+'motiflength.txt', 'w')
	f4.write(str(motif_length))
	f4.close();
	
		
		
			

	
main()
