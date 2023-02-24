import sys

def open_file1(path1,path2):

	dict_geneIDs = {}

	with open(path1, 'r') as f1:

		for line in f1.readlines():

			line = line.strip().split("\t")

			geneID = line[0].replace("ID=","")

			dict_geneIDs[geneID] = 1


	with open(path2, 'r') as f2:

		for line in f2.readlines():

			line = line.strip().split(" ")

			for geneID_ortho in line:

				if geneID_ortho in dict_geneIDs:
					print(str(geneID_ortho+"\t"+str(line).replace("[","").replace("]","").replace("'","").replace(",","")).replace(" ","\t"))


file1 = sys.argv[1]

open_file1(file1,"Orthogroups.txt")