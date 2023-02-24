import sys

def open_gtf_file(path,chrom):

	dict_genes = {}

	with open(path, 'r') as f1:

		for line in f1.readlines():

			line = line.strip().split("\t")

			key_id = line[8].split(";")

			if line[0] == chrom:

				dict_genes[key_id[0]] = [line[0],int(line[3]),int(line[4]),line[6],str(key_id[1:])]

	return dict_genes

def open_gff(path,chrom):

	dict_introns = {}

	with open(path, 'r') as f1:

		for line in f1.readlines():

			line = line.strip().split("\t")

			if line[5] not in dict_introns and line[0] == chrom:
				
				dict_introns[line[5]] = [[line[0],line[3],line[4]]]
			
			if line[5] in dict_introns and line[0] == chrom:
				
				dict_introns[line[5]].append([line[0],line[3],line[4]])				

	return dict_introns


def open_sv_file(path,dict_genes,dict_introns,chrom,svtype,trait):

	sep = "\t"

	out_sv_introns = open(svtype+"/outfile_"+chrom+"_"+svtype+"_"+trait+"_sv_introns_combinations.csv", 'w')
	out_sv_exons = open(svtype+"/outfile_"+chrom+"_"+svtype+"_"+trait+"_sv_exons_combinations.csv", 'w')
	out_5_utr = open(svtype+"/outfile_"+chrom+"_"+svtype+"_"+trait+"_sv_5utr_combinations.csv", 'w')
	out_3_utr = open(svtype+"/outfile_"+chrom+"_"+svtype+"_"+trait+"_sv_3utr_combinations.csv", 'w')
	out_intergenic = open(svtype+"/outfile_"+chrom+"_"+svtype+"_"+trait+"_sv_intergenic_control.csv", 'w')

	with open(path, 'r') as f1:

		for line in f1.readlines():

			line = line.strip().split("\t")

			#if line[0] != "Chrom":
			if  line[0] == chrom:

				sv_chrom = line[0]
				sv_start = int(line[1])

				switch_inter = 0
				for gene in dict_genes.keys():

					#Erster Fall: SV complete genic
					if sv_chrom == dict_genes[gene][0] and sv_start >= dict_genes[gene][1] and sv_start <= dict_genes[gene][2]:

						switch_inter = 1
						switch = 0

						if gene in dict_introns: #check if gene has introns

							for intron in dict_introns[gene]:

								if sv_chrom == intron[0] and sv_start >= int(intron[1]) and sv_start <= int(intron[2]):

									out_sv_introns.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
									switch = 1
									break

							if switch == 0:
								out_sv_exons.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
								break

					####5kb 3,5 UTRs

					#SV in 5kb start --> 5/3 UTR
					if sv_chrom == dict_genes[gene][0] and sv_start >= dict_genes[gene][1]-5000 and sv_start < dict_genes[gene][1]:

						switch_inter = 1
						if dict_genes[gene][3] == "+":
							out_5_utr.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
							break
						else:
							out_3_utr.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
							break
						
					#SV in 5kb stop --> 3/5 UTR
					if 	sv_chrom == dict_genes[gene][0] and sv_start > dict_genes[gene][2] and sv_start <= dict_genes[gene][2]+5000:

						switch_inter = 1					
						if dict_genes[gene][3] == "+":
							out_3_utr.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
							break
						else:
							out_5_utr.write(gene+"\t"+sep.join(line)+"\t"+dict_genes[gene][4]+"\n")
							break

				if switch_inter == 0:
					out_intergenic.write(sep.join(line)+"\t"+dict_genes[gene][4]+"\n")

	out_sv_introns.close()
	out_sv_exons.close()
	out_5_utr.close()
	out_3_utr.close()
	out_intergenic.close()

chrom = sys.argv[2] #chromosome subset
svtype = sys.argv[3] #SNPs
trait = sys.argv[4] #trait

dict_genes = open_gtf_file("Agria_Gene_Annotation_New/Agria_PASA_AGAT.all.genesRow.gff3",chrom)

dict_introns = open_gff("Agria_Gene_Annotation_New/Agria_PASA_AGAT.all.intronsRow.gff3",chrom)

open_sv_file(sys.argv[1],dict_genes,dict_introns,chrom,svtype,trait)
