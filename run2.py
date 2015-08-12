import sys
import subprocess
import random

def createpwm(i):
	absolute = 1
	padding = .01
	equivalent = .25
	Astring = str(absolute)+'\t'+str(padding)+'\t'+str(padding)+'\t'+str(padding)+'\n'
	Cstring = str(padding)+'\t'+str(absolute)+'\t'+str(padding)+'\t'+str(padding)+'\n'
	Tstring = str(padding)+'\t'+str(padding)+'\t'+str(absolute)+'\t'+str(padding)+'\n'
	Gstring = str(padding)+'\t'+str(padding)+'\t'+str(padding)+'\t'+str(absolute)+'\n'
	dString = str(equivalent)+'\t'+str(equivalent)+'\t'+str(equivalent)+'\t'+str(equivalent)+'\n'
	path = 'dataset'+str(i)+'/motif.txt'
	motif_file = open(path)
	data = motif_file.readline().split('\t')
	motif_file.close()
	outpath = 'dataset'+str(i)+'/pwm.txt'
	outfile = open(outpath, 'w')
	outfile.write('>'+data[0]+' '+data[1]+'\n')
	for row in range(int(data[1])):
		if data[2][row] is 'A':
			outfile.write(Astring)
		if data[2][row] is 'C':
			outfile.write(Cstring)
		if data[2][row] is 'T':
			outfile.write(Tstring)
		if data[2][row] is 'G':
			outfile.write(Gstring)
		if data[2][row] is '*':
			outfile.write(dString)
	outfile.write('<\n')
	outfile.close()
	
	
def main():
	total = 1
	num_datasets = 70
	if len(sys.argv) < 2:
		start = 1
	else:
		start = int(sys.argv[1])
	for i in range(start,num_datasets+1):
		createpwm(i);
		subprocess.call(['python', 'findmotif.py',str(i)])
main()
