from Bio import SeqIO
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq

ha_path_selected = "./H3N2/trimmed_msa_ha.fa"

path = "./../H3N2/removed_duplicates/"

first_bin = "06_07"
last_bin = "19_20"


new_file = "H3N2/labels_dna.fa"


def find_mutations(seq1, seq2, protein = "ha"):
	return [protein + "_" + seq1[i - 1] + str(i) + seq2[i - 1] for i in range(1, len(seq1) + 1) if seq1[i - 1] != seq2[i - 1] and seq1[i - 1] != "-" and seq2[i - 1] != "-"]


def padding_zero(i):
	return str(i).zfill(2)
	
list_key = [seq_record.id for seq_record in SeqIO.parse(ha_path_selected, "fasta")]
list_value = [element for nestedlist in [15 * ["{}_{}".format(padding_zero(i), padding_zero(i+1))] for  i in range(int(first_bin[0:2]), int(last_bin[0:2]) + 1)] for element in nestedlist]


sequence_name_dict = dict(zip(list_key, list_value))

key_flu_season = {}
for name, flu_season in sequence_name_dict.items():
	if flu_season not in key_flu_season:
		key_flu_season[flu_season] = [name]
	else:
		key_flu_season[flu_season].append(name)


all_proteins_sequences = {}
for flu_season in key_flu_season:
	path_ha_proteins = "{}/ha_proteins/ha_{}.fa".format(path, flu_season)
	path_na_proteins = "{}/na_proteins/na_{}.fa".format(path, flu_season)
	for seq_record in SeqIO.parse(path_ha_proteins, "fasta"):
		if seq_record.id in key_flu_season[flu_season]:
			all_proteins_sequences[seq_record.id] = {}
			all_proteins_sequences[seq_record.id]["ha_sequence"] = seq_record
	for seq_record in SeqIO.parse(path_na_proteins, "fasta"):
		if seq_record.id in key_flu_season[flu_season]:
			all_proteins_sequences[seq_record.id]["na_sequence"] = seq_record

list_of_positions_ha = [158, 213, 277]
list_of_positions_na = []


labels_dict = {}
for seq_record in all_proteins_sequences:
	
	ha_seq = all_proteins_sequences[seq_record]["ha_sequence"]
	na_seq = all_proteins_sequences[seq_record]["na_sequence"]
	aminoacid_position_ha = "".join([ha_seq[position - 1] for position in list_of_positions_ha])
	labels_dict[seq_record] = aminoacid_position_ha
	aminoacid_position_na = "".join([na_seq[position - 1] for position in list_of_positions_na])
	labels_dict[seq_record] += aminoacid_position_na

labels = []
for seq_record in all_proteins_sequences:
	
	all_proteins_sequences[seq_record]["ha_sequence"].seq = Seq(labels_dict[seq_record])
	labels.append(all_proteins_sequences[seq_record]["ha_sequence"])

SeqIO.write(labels, new_file, "fasta")