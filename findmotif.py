import numpy as np
from numpy import *
import random
import os
from scipy import stats
import sys
import time

start_time = time.time()

def lists_equal(l1, l2):
	for i in range(len(l1)):
		if l1[i] != l2[i]:
			return False
	return True

# get motif
setnum = sys.argv[1]
path = "dataset"+setnum+"/"
print path
motif_length_file = open(path+'motiflength.txt','r')
motif_length = int(motif_length_file.read())
print('motif length = '+str(motif_length))
motif_length_file.close()
# get sequences
sequences = []
sequence_length = 0

seq_file = open(path+'sequences.fa','r')
for line in seq_file:
	if ">" not in line:
		sequences.append(line)

seq_file.close()

num_sequences = len(sequences)
sequence_length = len(sequences[0])-1
print('Number of sequences = '+str(num_sequences)+". Sequence length = "+str(sequence_length))

outer_iterations = 500
best_start_positions = np.zeros((outer_iterations, 1+num_sequences)) # extra column for the score of this set

for k in range(outer_iterations):
	print('dataset '+setnum+', start positions '+str(k))
	not_improved = 0  	# when this reaches 50, break out of the while loop b/c it means the score is unlikely to improve anymore
	old_score = 0
	new_score = 0

	# initialize random starting positions:
	start_positions = []
	for i in range(num_sequences):
		start_positions.append(random.randint(0, sequence_length-motif_length))

	while True:
		# Step 1: Get background frequencies (https://www.biostat.wisc.edu/bmi776/spring-07/lectures/gibbs-example.pdf)
		background_frequencies = np.zeros((4,))
		for i in range(num_sequences):
			seq = sequences[i]
			for j in range(sequence_length):
				if j < start_positions[i] or j >= start_positions[i]+motif_length:
					if seq[j] == 'A':
						background_frequencies[0] += 1
					if seq[j] == 'C':
						background_frequencies[1] += 1
					if seq[j] == 'T':
						background_frequencies[2] += 1
					if seq[j] == 'G':
						background_frequencies[3] += 1

		background_frequencies /= num_sequences*(sequence_length-motif_length)

		# Step 2: Choose one of the sequences at random to remove
		removed = random.randint(0,num_sequences-1)
		#print 'removed sequence {}'.format(removed)
		# get the lmers starting at the random positions in the remaining sequences
		substrings = []
		for i in range(num_sequences):
			if i != removed:
				start = start_positions[i]
				seq = sequences[i]
				substring = seq[start:start+motif_length]
				substrings.append(substring)
				#print 'substring[{}] \t {}'.format(i,substring)

		# Step 3: Create profile P from l-mers in remaining sequences
		profile_matrix = np.zeros((4, motif_length))
		for col in range(motif_length):
			for substring in substrings:
				letter = substring[col]
				if letter == 'A':
					profile_matrix[0,col] += 1
				elif letter == 'C':
					profile_matrix[1,col] += 1
				elif letter == 'T':
					profile_matrix[2,col] += 1
				elif letter == 'G':
					profile_matrix[3,col] += 1		
			# divide every entry in that column by the total of the column (= num_sequences-1):
			profile_matrix[:,col] /= (num_sequences-1)
		#print profile_matrix

		# Step 4: Calculate Pr(a|P) for every possible 8-mer in the removed sequence
		lmer_probabilities = []
		non_zero_positions = []
		lowest_probability = 1 # need this for Step 5
		removed_sequence = sequences[removed]

		for i in range(sequence_length-motif_length+1):
			lmer = removed_sequence[i:i+motif_length]
			probability = 1
			for j in range(motif_length):	
				# for each letter in the l-mer, retrieve its probability from the profile matrix and update probability
				c = lmer[j]
				if c == 'A':
					probability *= profile_matrix[0,j]
				elif c == 'C':
					probability *= profile_matrix[1,j]
				elif c == 'T':
					probability *= profile_matrix[2,j]
				elif c == 'G':
					probability *= profile_matrix[3,j]

			if probability != 0:
				#print probability
				lmer_probabilities.append(probability)
				non_zero_positions.append(i)
				if probability < lowest_probability:
					lowest_probability = probability

		# Step 5: divide each probability Pr(a|P) by the lowest probability to get the whole-number ratio
		denominator = 0
		ratio = list(lmer_probabilities)
		for i in range(len(lmer_probabilities)):
			ratio[i] /= lowest_probability
			denominator += ratio[i]
		
		# Define probabilities of starting positions according to the computed ratios. 
		xk = non_zero_positions
		pk = [x / denominator for x in ratio]
		if len(pk) != 0:
			custm = stats.rv_discrete(name='custm', values=(xk,pk))
			potential_start = custm.rvs()

		# Recalculate score of the lmer in the removed sequence
		new_score = 0
		removed_lmer = sequences[removed][potential_start:potential_start+motif_length]
		for i in range(motif_length):
			c = removed_lmer[i]
			r = 0 
			if c == 'A':
				r = 0
			elif c == 'C':
				r = 1
			elif c == 'T':
				r = 2
			elif c == 'G':
				r = 3		
			new_score += (profile_matrix[r,i]/background_frequencies[r])

		# if new score >= old_score, jump to that new start position.
		# otherwise, keep the old position
		if new_score >= old_score:
			if new_score - old_score == 0:
				not_improved += 1
			else:
				not_improved = 0
			start_positions[removed] = potential_start
			old_score = new_score
		else:
			not_improved += 1
		if not_improved == 50:
			break

	# save this set of start positions and its score in best_start_positions:
	best_start_positions[k][0] = old_score

	for c in range(1, len(best_start_positions[k])):
		best_start_positions[k][c] = start_positions[c-1]

#print best_start_positions
best_score_idx = 0
for i in range(outer_iterations):
	if best_start_positions[i][0] > best_start_positions[best_score_idx][0]:
		best_score_idx = i

final_start_positions = best_start_positions[best_score_idx][1:]

elapsed_time = time.time() - start_time
time_file = open(path+'time.txt', 'w')
time_file.write('elapsed time: '+str(elapsed_time)+' seconds')
time_file.close()

sites_file = open(path+'predictedsites.txt','w')
motif_file = open(path+'predictedmotif.txt','w')
lmers_file = open(path+'lmers.txt','w')
predicted_motif = np.zeros((motif_length, 4))
lmers = []

for i in range(num_sequences):
	seq = sequences[i]
	start = (int)(final_start_positions[i])
	sites_file.write('S'+str(i)+'\t'+str(start)+'\n')
	lmer = seq[start:start+motif_length]
	print(lmer)
	lmers.append(lmer)
	for j in range(motif_length):
		c = lmer[j]
		if c == 'A':
			predicted_motif[j][0] += 1
		elif c == 'C':
			predicted_motif[j][1] += 1
		elif c == 'T':
			predicted_motif[j][2] += 1
		elif c == 'G':
			predicted_motif[j][3] += 1

lmers_file.write("\n".join(lmers))
motif_file.write('>PMOTIF'+'\t'+str(motif_length)+'\n')
for row in range(motif_length):
	line = predicted_motif[row]
	s = str(int(line[0]))+'\t'+str(int(line[1]))+'\t'+str(int(line[2]))+'\t'+str(int(line[3]))+'\n'
	motif_file.write(s)

motif_file.write('<\n')
lmers_file.close()
sites_file.close()
motif_file.close()
print('dataset'+ str(setnum)+' is done')
print('-------------------------------------------------')
