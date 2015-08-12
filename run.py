#!/usr/bin/python2.7
import sys
import subprocess
import random

def main():
	total = 1
	num_datasets = 70
	for i in range(num_datasets):
		NM=1		# number of variable positions
		ML=8		# motif length
		SL=500		# sequence length
		SC=10		# sequence count
		if i >= 10 and i < 20:
			NM=0
		elif i >= 20 and i < 30:
			NM=2
		elif i >= 30 and i < 40:
			ML=6
		elif i >= 40 and i < 50:
			ML=7
		elif i >= 50 and i < 60:
			SC=5
		elif i >= 60 and i <= 70:
			SC=20
		print str(ML)+" "+str(NM)+" "+str(SL)+" "+str(SC)+" "+str(total)
		subprocess.call(['python', 'part1.py',str(ML),str(NM),str(SL),str(SC),str(total)])
		#subprocess.call(['python', 'findmotif.py',str(i)])	# pass into findmotif.py the dataset #
		total=total+1



main()
