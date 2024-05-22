## 2021 05 24
## Dr David Wright
## Parse GTEx excel sheet of individuals/tissues/SRR accessions and the splice counts file to calculate PSI for tissue-specific exon
## source python-3.5.1
import sys
import csv 
import pandas as pd
import numpy 
import scipy.stats

# inputs
gtex_file = sys.argv[1] # count file
searchlist = sys.argv[2] # SRA list
tissue_options = sys.argv[3] # comma sep text file of tissue options to narrow the GTEx to e.g. 'brain,aorta'
gene = sys.argv[4] # gene identifiers that you want to extract such as for kalirin eg. 'ENSG00000160145'
tissueA = sys.argv[5] # text for comparison tissue1 i.e. "Brain" to be used in all produced results 
tissueB = sys.argv[6] # string of tissue2 e.g. Aorta, or Oesophagus to be used in all produced results 

def main(searchlist, tissue_options):
	with open(searchlist, 'r') as sra:
		sradict = {}
		for line in sra: 
			array = line.split('\t')
			#print(array[24], array[25], array[20], array[26]) # 20 is the SRR, 24 is the GTEX 16 char and 25 is the tissue 39 is subjectID
			if array[24] not in sradict:
				sradict[array[24]] = tuple([array[25],array[39]])
		# read in the tissue options to filter the dict by
		# print(sradict)
		tissue_file = open(tissue_options, 'r')
		tissues = tissue_file.readline().strip().split(',')
		print(tissues)
		gtex_tissue = []
		gtex_key = []
		outfile = open('GTEx_accessions_for_splicechecks_'+str(tissueA)+'_'+str(tissueB)+'.txt', 'w')
		for k,v in sradict.items(): 
			for j in v:
				if any(substring in j for substring in tissues):
					gtex_tissue.append(k.strip()) # make list of GTEx identifiers to filter target file
					outfile.write(k+'\t'+str(v[0])+'\t'+str(v[1])+'\n') #can write this out if a GTExID --> tissue key is needed
					gtex_key.append([k,v])
	#print(gtex_tissue)
	print(str(len(gtex_tissue))+' columns matching requested tissue(s) by substring match')
	#print(gtex_tissue)
	outfile.close()
	return gtex_tissue, gtex_key

def split_gtex(gtex_file, gene, gtex_tissue):
	# read in huge file and split down to gene of interest!
	with open(gtex_file, 'r') as gtex:
		#next(gtex) # version
		#next(gtex) # row/column info
		header = next(gtex) # stick column headers in holder
		header_row = numpy.array(header.strip().split('\t')) # convert list to array for DF 
		#print(header_row)
		keep = []
		#keep.append(header.rstrip().split('\t'))
		for line in gtex: 
			if str(gene) in line.strip().split('\t')[1]:  # subset file by gene of choice, note .version number hence search 'in' & not 'equal to'...
				keep.append(line.strip().split('\t'))
	#print(keep)
	print('there are '+str(len(keep))+' junction entries in total for requested gene(s): '+str(gene))	
	gtexDF = pd.DataFrame(keep, columns = header_row)
	print('pre-filtering df size is '+str(gtexDF.shape))
	print(gtexDF.head())
	# FILTER & MERGE IDs AND FILTERED COUNTS: filter the gtexDF to only include the columns of interest 
	# split counts and IDs to filter counts by list of matching GTEx codes (which would remove ID columns!)
	gtexDF_counts = gtexDF.drop(gtexDF.columns[[0,1]], axis=1)
	gtexDF_IDs = gtexDF.filter(gtexDF.columns[[0,1]], axis=1)
	#Filter by the list, using list comprehension to avoid issues with missing keys 
	gtexDF_filt = gtexDF_counts.drop(columns=[col for col in gtexDF_counts if col not in gtex_tissue])
	#Merge first columns back to counts 
	gtexsubset =  pd.concat([gtexDF_IDs, gtexDF_filt], axis=1)
	print('after filtering the df is '+str(gtexsubset.shape))
	print(gtexsubset.head())
	gtexsubset.to_csv('subsetted_junctions_'+str(tissueA)+'_'+str(tissueB)+'.txt', index=False)
	return gtexsubset

def calculate_PSI(tissueA, tissueB, gtex_counts, gtex_key):
	#identify the junction of interest & tissues of interest
	# find pairs of samples to use in the comparison 
	comp_dict = {}
	for line in gtex_key:
		(tissue, subject) = line[1]
		if not subject in comp_dict: 
			comp_dict[subject] = [tuple([line[0], tissue])]
		else:
			comp_dict[subject].append(tuple([line[0], tissue]))
	# find paired samples where each tissue is represented
	subjects = []
	tissue1 = []
	tissue2 = []
	for k,v in comp_dict.items(): #this asserts only multiple datasets per individual 
		if len(v) >= 2:
			if any(str(tissueB) in element for entry in v for element in entry) and any(str(tissueA) in element for entry in v for element in entry): #this asserts at least one of each tissue type per individual
				#print(k,v)
				subjects.append(k)
	# extract the sample IDs that need the counts for these individuals, now that they are paired (irrespective of how many of each tissue is there)
	for k,v in comp_dict.items():
		if k in subjects:
			for element in v:
				if str(tissueB) in element[1]:
					tissue2.append(tuple([k, element[0]]))
				elif str(tissueA) in element[1]:
					tissue1.append(tuple([k, element[0]]))
	#print('the list of aorta accessions: '+ str(tissue2)+'\n')
	#print('the list of brain accessions: '+str(brain)+'\n')
	#print(subjects)
	#print(tissue1) #brain
	#print(tissue2) # other eg aorta
	# extract relevant counts for each set of reads. For this the exon details are: 
	
	# KALRN EXON DETAILS
	#Exon before: chr3:124,633,849-124,633,953;
	#Brain-specific exon: chr3:124,637,208-124,637,303;
	#Exon after: chr3:124,650,808-124,650,938
	# So the STAR intron-based counts mean the corresponding reads needed are: 
	# inclusion: A) chr3_124633954_124637207 and B) chr3_124637304_124650807
	# exclusion: C) chr3_124633954_124650807
	incl_juncs = ['chr3_124633954_124637207', 'chr3_124637304_124650807']
	excl_juncs = 'chr3_124633954_124650807'
	tissue1_incl = []
	tissue1_excl = []
	tissue2_incl = []
	tissue2_excl = []
	ind_counts = {} #individual counts across the four groups
	# set Name column as index for searching the intron coordinates
	counts = gtex_counts.set_index('Name')
	#print(str(counts.index.tolist()))
	missingGTExfromSra = []
	check_cols = list(counts)
	
	## GET BRAIN COUNTS - AND MERGE INDIVIDUALS 
	for ID in tissue1:  # extracting the correct element of the tuple (this 16chr-GTEx ID rather than subject ID)
		print(ID)
		if ID[1] in check_cols:
			if ID[0] not in ind_counts.keys():
				ind_counts[ID[0]] = {str(tissueA)+'_incl': [], str(tissueA)+'_excl': []} #make entry for individual 
			
			for junc in incl_juncs:
				tissue1_incl.append(counts.loc[str(junc).strip(),str(ID[1])])
				ind_counts[ID[0]][str(tissueA)+'_incl'].append(int(counts.loc[str(junc).strip(),str(ID[1])]))
			tissue1_excl.append(counts.loc[str('chr3_124633954_124650807'),str(ID[1])])
			ind_counts[ID[0]][str(tissueA)+'_excl'].append(int(counts.loc[str('chr3_124633954_124650807'),str(ID[1])]))
		elif ID[1] not in check_cols:
			missingGTExfromSra.append(ID[1]) 
	#print(ind_counts)	
	## GET SECOND TISSUE COUNTS -- AND MERGE INDIVIDUALS. ALL KEYS SHOULD EXIST AS PAIRED SAMPLES!!!
	for ID in tissue2:
		if ID[1] in check_cols: 
			if ID[0] not in ind_counts.keys():
				print('second tissue has data/sample '+str(ID[0])+' '+str(ID[1])+' which is not found in brain set meaning an sra file - GTEx count file conflict. Check "Sra_entries_missing_in_GTExcount_file.txt" for culprit!')# enforce paired-only by removing sra/GTEx list conflicts (ie. samples removed from brain as not present in GTEx despite presence in SRA...then remove corresponding second tissue's entries that match that ID!)
			else:
				ind_counts[ID[0]][str(tissueB)+'_incl'] = []
				ind_counts[ID[0]][str(tissueB)+'_excl'] = [] 
			if ID[0] in ind_counts.keys(): 
				for junc in incl_juncs:
					tissue2_incl.append(counts.loc[str(junc).strip(),str(ID[1])])
					ind_counts[ID[0]][str(tissueB)+'_incl'].append(int(counts.loc[str(junc).strip(),str(ID[1])]))
				tissue2_excl.append(counts.loc[str('chr3_124633954_124650807'),str(ID[1])])
				ind_counts[ID[0]][str(tissueB)+'_excl'].append(int(counts.loc[str('chr3_124633954_124650807'),str(ID[1])]))
		elif ID[1] not in check_cols:
			missingGTExfromSra.append(ID[1])
	
	print(ind_counts)
	# conver to int for summing!
	tissue1_incl = [int(i) for i in tissue1_incl]
	tissue1_excl = [int(i) for i in tissue1_excl]
	tissue2_incl = [int(i) for i in tissue2_incl]
	tissue2_excl = [int(i) for i in tissue2_excl]
	
	BI = sum(tissue1_incl)
	print(str(tissueA)+' incl. exon: '+ str(BI))
	BE = sum(tissue1_excl)
	print(str(tissueA)+' excl. exon: '+ str(BE))
	BTOT = BI+BE
	print(str(tissueA)+' total reads: '+ str(BTOT))
	TI = sum(tissue2_incl)
	print(str(tissueB)+' incl. exon: '+ str(TI))
	TE = sum(tissue2_excl)
	print(str(tissueB)+' excl. exon: '+ str(TE))
	TTOT = TI+TE
	print(str(tissueB)+' total reads: '+ str(TTOT))
	
	### WRITE OUT SRA-->GTEX CONFLICTS
	missing = open('Sra_entries_missing_in_GTExcount_file_'+str(tissueA)+'_'+str(tissueB)+'.txt', 'w')
	for element in missingGTExfromSra:
		missing.write(str(element)+'\n') #write out the mismatch list between Sra and GTEx file for reference...

	## FISHER EXACT TEST of exon inclusion, normalising for total reads:
	oddsratio, pval = scipy.stats.fisher_exact([[BI,(BTOT-BI)],[TI,(TTOT-TI)]])
	print('odds ratio: '+str(oddsratio))
	print('p-value: '+str(pval))
	
	### PRINT OUT THE INDIVIDUAL COUNTS FROM ind_counts DICT
	# sum up/merge the counts for each subjects' tissue category
	for key,val in ind_counts.items():
		for tissue,counts in val.items():
			val[tissue] = sum(counts)
	fields = [ 'subjectID',str(tissueA)+'_incl',str(tissueA)+'_excl',str(tissueB)+'_incl', str(tissueB)+'_excl']
	with open('subject_counts_'+str(tissueA)+'_'+str(tissueB)+'.csv', 'w') as f:
		w = csv.DictWriter(f,fields)
		w.writeheader()
		for key,val in sorted(ind_counts.items()):
			row = {'subjectID': key}
			row.update(val)
			w.writerow(row)
		
# RUN PROGRAM
if __name__ == '__main__':
	# read in the input vars and produce list of GTEX codes that the file needs slicing to
	gtex_tissue,gtex_key = main(searchlist, tissue_options)
	#print(gtex_key)
	gtex_counts = split_gtex(gtex_file, gene, gtex_tissue)
	calculate_PSI(tissueA, tissueB, gtex_counts, gtex_key)


