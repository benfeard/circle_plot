'''
by Benfeard Williams
Web Feb 8 2017

This script processes the raw input for cross correlation calculations from Wordom.

usage: process_raw_crosscorr.py crosscorr.wordom outputFilename
'''
import sys

output = open(sys.argv[-1][:-3] + "final.dat", 'w')

with open(sys.argv[-1]) as myFile:
	for line in myFile:
		if '#' not in line:
			line = line.strip().split()
			res1 = int(line[0])
			res2 = int(line[1])
			#if abs(res1-res2) <= 15: continue
			corr = float(line[-1])
			output.write("%d\t%d\t%f\n" % (res1,res2,corr))