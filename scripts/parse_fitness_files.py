#!/usr/bin/env python
#notes for parsing the viral fitness files

from __future__ import print_function
import argparse
import os
from datetime import date
import re
from pprint import pprint as pp
import csv
import sys



parser = argparse.ArgumentParser(description='''
The point of this script is to get the previously published fitness data into a usable format.

This script takes a text file and parses it into columns

By Shawn Whitefield 2-2-2016

 ''')

parser.add_argument('-f', help='input file')
parser.add_argument('-o', help = 'output file')
args=parser.parse_args()

# open the files
infile = open(args.f, 'r') # file to be parsed
outfile = open(args.o, 'w') # output file

for line in infile:
		fields= re.sub('[()]',"",line)
		fields = re.sub('lethal',"NA NA",fields)
		test = re.split(r'(.*?\s.*?\s.*?)\s', fields)
		test1 = str(test)
		test = test1.replace("[" , "")
		test = test.replace("]","")
		test = test.replace(",","\n")
		test = test.replace("'",'')
		test = test.replace(" ", '\t')
		print(test)

		
		#print(test)
		#test = '\n'.join(test1)
		#print(test)
		#test = test.replace('. '',',test)
		#print(test)









