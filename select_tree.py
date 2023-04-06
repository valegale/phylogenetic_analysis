from Bio import SeqIO
import random
import csv

###  parameters 

bin_start = "06_07"
bin_end = "19_20"
n_sequences = 15
flu_season_outlier = "04_05"
year_outlier = "2004"

def select_random_proteins(path_sequences, N, first_bin, last_bin):
	
	def padding_zero(i):
		return str(i).zfill(2)
	
	total_list_sequences_ha = []
	total_list_sequences_na = []
	name_year = []
	for i in range(int(first_bin[0:2]), int(last_bin[0:2]) + 1):
		bin_number = "{}_{}".format(padding_zero(i), padding_zero(i+1))
		
		
		ha_protein_file_bin = "{}/ha_proteins/ha_{}.fa".format(path_sequences, bin_number)
		na_protein_file_bin = "{}/na_proteins/na_{}.fa".format(path_sequences, bin_number) 
		
		sequences_bin = {}
		for seq_record in SeqIO.parse(ha_protein_file_bin, "fasta"):
			sequences_bin[seq_record.id] = {}
			sequences_bin[seq_record.id]["ha"] = seq_record
			
		for seq_record in SeqIO.parse(na_protein_file_bin, "fasta"):
			sequences_bin[seq_record.id]["na"] = seq_record
		
		sequences_id = random.sample(list(sequences_bin.keys()), N)
		
		for sequence_id in sequences_id:
			total_list_sequences_ha.append(sequences_bin[sequence_id]["ha"])
			total_list_sequences_na.append(sequences_bin[sequence_id]["na"])
			name_year.append((sequence_id, bin_number))

	return total_list_sequences_ha, total_list_sequences_na, name_year

def adding_outlier(selected_sequences_ha, selected_sequences_na, outlier_ha, outlier_na):
	for seq_record in SeqIO.parse(outlier_ha, "fasta"):
		name_outlier = seq_record.id
		selected_sequences_ha.append(seq_record)
	for seq_record in SeqIO.parse(outlier_na, "fasta"):
		selected_sequences_na.append(seq_record)
	return selected_sequences_ha, selected_sequences_na, name_outlier
	
	
def clean_ambiguous_characters(selected_sequences_ha, selected_sequences_na):

	for sequence in selected_sequences_ha:
		for i in range(len(sequence)):
			if sequence[i] == "J":
				print (i)
				sequence.seq = sequence.seq[:i] + "X" + sequence.seq[i + 1:]
	for sequence in selected_sequences_na:
		for i in range(len(sequence)):
			if sequence[i] == "J":
				sequence.seq = sequence.seq[:i] + "X" + sequence.seq[i + 1:]	
	return selected_sequences_ha, selected_sequences_na



def extract_nucleotide_data(path_sequences, name_year):
	total_list_sequences_ha = []
	total_list_sequences_na = []
	
	sequences_bin = {}
	for (name, flu_season) in name_year:
		
		ha_dna_file_bin = "{}/ha_filtered/ha_{}.fa".format(path_sequences, flu_season)
		na_dna_file_bin = "{}/na_filtered/na_{}.fa".format(path_sequences, flu_season) 
		
		for seq_record in SeqIO.parse(ha_dna_file_bin, "fasta"):
			if seq_record.id == name:
				sequences_bin[seq_record.id] = {}
				sequences_bin[seq_record.id]["ha"] = seq_record
			
		for seq_record in SeqIO.parse(na_dna_file_bin, "fasta"):
			if seq_record.id == name:
				sequences_bin[seq_record.id]["na"] = seq_record
		
	for sequence_id in sequences_bin:
		total_list_sequences_ha.append(sequences_bin[sequence_id]["ha"])
		total_list_sequences_na.append(sequences_bin[sequence_id]["na"])
	print (len(sequences_bin))
	return total_list_sequences_ha, total_list_sequences_na

path_sequences = "./../H3N2/removed_duplicates"	
path_sequences_nucleotides = "./../H3N2"
path_results = "H3N2"
outlier_ha = "{}/outlier_ha.fa".format(path_results)
outlier_na = "{}/outlier_na.fa".format(path_results)


selected_sequences_ha, selected_sequences_na, name_year = select_random_proteins(path_sequences, n_sequences, bin_start, bin_end)

#print (selected_sequences_ha)
#optional, create a nucleotide file
dna_sequences_ha,dna_sequences_na = extract_nucleotide_data(path_sequences_nucleotides, name_year)



selected_sequences_ha, selected_sequences_na = clean_ambiguous_characters(selected_sequences_ha, selected_sequences_na)
selected_sequences_ha, selected_sequences_na, name_outlier = adding_outlier(selected_sequences_ha, selected_sequences_na, outlier_ha, outlier_na)


# silenced to avoid casual running and overwriting of the sequences in the folder
#SeqIO.write(selected_sequences_ha, "{}/selected_sequences_ha.fa".format(path_results), "fasta")
#SeqIO.write(selected_sequences_na, "{}/selected_sequences_na.fa".format(path_results), "fasta")
# dna data
#SeqIO.write(dna_sequences_ha, "{}/selected_sequences_dna_ha.fa".format(path_results), "fasta")
#SeqIO.write(dna_sequences_na, "{}/selected_sequences_dna_na.fa".format(path_results), "fasta")


'''
with open('{}/metadata.tsv'.format(path_results), 'w', newline='') as tsvfile:
	writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
	writer.writerow(["sequence_id", "flu_season", "year"])
	writer.writerow([name_outlier, flu_season_outlier, year_outlier])
	for (name, flu_season) in name_year:
		writer.writerow([name, flu_season, name[-4:]])
	'''
	
	