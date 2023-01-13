import sys
import os


list_files = os.listdir()

for file in list_files:

	if file[0:8]  == "QTLs_sig":

		dict_add = {}
		dict_gen = {}

		with open(file, 'r') as f1:
		
			for line in f1.readlines():
				
				line = line.strip().split(",")

				if len(line) > 1:

					if "additive" in line[2]:

						dict_add[line[4]] = 1

					else:
					
						dict_gen[line[4]] = 1

		count = 0

		for key in dict_add:

			if key in dict_gen:

				count+=1

		print(file,count)
