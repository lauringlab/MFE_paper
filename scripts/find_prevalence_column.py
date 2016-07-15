from __future__ import print_function
import re
import sys



#THIS FINDS THE LINES THAT MATCH THE AMINO ACID SUB THAT WE HAVE IN OUR COLLECTION
#PRINTS A COLUMN THAT CAN BE CBINDED TO R DATAFRAME 
#DONT INCLUDE HEADERLINE (FIRST NA)	
	
with open(sys.argv[1]) as fp:
	with open(sys.argv[2],'w') as out:
		next(fp) #skip the header
		for line in fp:
			line = line.strip()
			line = line.split('\t')
			str = line[8]
			search_string = str + "=.*?,"
			match = re.search(search_string, line[2])
			if match:
				print(match.group(0))
				out.write(match.group(0)+"\n")
			else:
				print("NA")
				out.write("NA"+"\n")
