import sys
import collections

def open_variants_file(path,distance_variants):

	dict_common_variants = {}

	with open(path, 'r') as f1:

		help_list_start = []
		c = 1 #window counter

		for line in f1.readlines():
			line = line.strip().split("\t")

			#line[1] = start position
			if dict_common_variants == {}: #to avoid an error message for the first item

				#should only jump for the first variant into this if statement
				dict_common_variants[1] = [[int(line[1]),int(line[2]),line[3]+"_"+line[4]]] 
				help_list_start.append(int(line[1]))

			else:

				help_list_start.append(int(line[1])) #safe all start positions in this list
				help_counter = 0

				for i in range(0,distance_variants):

					if int(line[1])+i == help_list_start[-2] or abs(help_list_start[-2] - (int(line[1])+i)) <= distance_variants: #similar start position or difference smaller than 50bp

							
						dict_common_variants[c].append([int(line[1]),int(line[2]),line[3]+"_"+line[4]])
						break

					else:
						help_counter += 1 #count non common start positions

				if help_counter == distance_variants:

					c += 1
					dict_common_variants[c] = [[int(line[1]),int(line[2]),line[3]+"_"+line[4]]]

	return dict_common_variants #Hash includes windows with similar and  <50bp start positions to the start position of the next variant
def control_stopInCluster(dict_starts,distance_variants):

	dict_newClusters = {} #Safe clusters which should be removed because of her stop position
	cw = 1 #new window counter for this dict

	for clusters in dict_starts.keys():

		if len(dict_starts[clusters]) > 1:

			clustersSortedEnd = sorted(dict_starts[clusters],key=lambda x: x[1]) #sorts each list of lists by stop position

			help_list_stops = [] #stop position for each cluster 

			list_index_counter = 0 #count indices of clusters to remove from help_list_stop

			checkNewCluster = [] #cluster which will be safed in new dict

			for cluster in clustersSortedEnd:

				if list_index_counter == 0:
					help_list_stops.append(cluster[1])
				else:
					help_list_stops.append(cluster[1])

					help_counter = 0

					for i in range(0,distance_variants):

						if cluster[1]+i == help_list_stops[-2] or abs(help_list_stops[-2] - (cluster[1]+i)) <= distance_variants: #similar stop position or difference smaller than 200bp

							#cluster bleibt hier erhalten
							break

						else:
							help_counter += 1 #count non common start positions

					if help_counter == distance_variants:

						#1) variant will be removed						
						dict_starts[clusters].remove(cluster)
					
						#2) variant will be removed from help list
						help_list_stops.pop(list_index_counter)
						list_index_counter -= 1
						
						#3) New entry in new dict
						
						checkNewCluster.append(cluster)
						cw += 1

				list_index_counter += 1

			if checkNewCluster != []:

				dict_newClusters[cw] = checkNewCluster #cw = cluster window


	return dict_starts, dict_newClusters
def genotyping(list_dict1,chromosome,path_outfile,list_samples):

	dict_allClusters = {}

	window = 1
	for run in list_dict1:
		for cluster in run.keys():
			dict_allClusters[window] = run[cluster]
			window += 1

	for cluster in dict_allClusters.keys():
		
		#alle ersten elemente sortieren und max and vorkommen + end
		help_start = []
		help_stop = []

		for start,stop,geno in dict_allClusters[cluster]:
			help_start.append(start)
			help_stop.append(stop)

		counter_start = collections.Counter(help_start)
		counter_stop = collections.Counter(help_stop)
		
		dict_allClusters[cluster].append([counter_start.most_common(1)[0][0]]) #as list
		dict_allClusters[cluster].append([counter_stop.most_common(1)[0][0]]) #as list ???

		#print(dict_allClusters[cluster])

			
		for sample_name in list_samples:

			if sample_name in str(dict_allClusters[cluster]):

				position_name = str(dict_allClusters[cluster]).find(sample_name)

				dosage = str(dict_allClusters[cluster])[position_name+len(sample_name)-1+2]

				dict_allClusters[cluster].append(dosage)
			else:
				dict_allClusters[cluster].append(0)
			
			##if sample_name not in str(dict_allClusters[cluster]):
			##	dict_allClusters[cluster].append(0)


	with open(path_outfile, 'w') as out:

		for final_clusters in dict_allClusters.keys():
			out.write(chromosome+"\t"+str(dict_allClusters[final_clusters][-110:]).replace("[","").replace("]","").replace("'","").replace(", ","\t")+"\n")

#input
distance_variants = int(sys.argv[2])
dict_starts = open_variants_file(sys.argv[1],distance_variants)
chromosome = sys.argv[3]
outfile = sys.argv[4]

#output
dict2 = "" #just to define for the first time
list_dict1 = []
list_dict2 = [dict_starts]

i = 0
while(dict2 != {}):

	dict1,dict2 = control_stopInCluster(list_dict2[i],distance_variants)
	list_dict1.append(dict1)
	list_dict2.append(dict2)
	i += 1

list_names = ["Adretta","Agria","Albtras","Allians","Altus","Ambition","Atlntic","Atzimba","belena","Bintje","BNA_1","BNA_2","BNA_3","BNA_4","BNA5","Cara","Celtane","CGN17881","CGN_18114","Charlte","Cherie","Colomba","Dark","Desiree","Donata","EaryRse","Edison","Eurogrande","Europrim","Felsina","Flava","Fontane","Gala","Gladitor","GLKS3","Granola","H98D12","Harpun","hermes","Innovat","Jelly","jukijiro","Karelia","karolin","Kathadin","KENN","KingRust","Kolibri","Krone","Kuba","Kuras","LadyRose","Laura","Leyla","Lilly","Marabel","MarsPipr","Natalia","Nevsky","Nicola","Odysseus","Olympus","ONA","Otolia","P3","P40","PentDell","Pirol","Princess","p_russet","Quadriga","Quarta","record","Regina","Rode_Est","Rooster","Rosara","Rudolph","russet_burbank","s13_017","S14214","S14317_1","S14_9122","S5061","S5089","S5202","S5220a","S5342","SA222_2","SA225_2","saskia","Semlo","seresta","shc909","Shepody","Skawa","snowden","solist","Spunta","Talent","tyoshi","Udacha","Velox","verdi","Vitablla","vr808","Yangana","zorba"]

genotyping(list_dict1,chromosome,outfile,list_names)

#Das sieht sehr gut aus!
#for run in list_dict1:
#	for cluster in run.keys():
#		print(cluster,run[cluster])
#######################################


