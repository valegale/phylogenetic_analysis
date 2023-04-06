from Bio import SeqIO
import pandas as pd
import seaborn as sns

ha_path_selected = "./H3N2/trimmed_msa_ha.fa"

path = "./../H3N2/removed_duplicates/"

first_bin = "06_07"
last_bin = "19_20"

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
		all_proteins_sequences[seq_record.id] = {}
		all_proteins_sequences[seq_record.id]["ha_sequence"] = seq_record.seq
	for seq_record in SeqIO.parse(path_na_proteins, "fasta"):
		all_proteins_sequences[seq_record.id]["na_sequence"] = seq_record.seq


flu_season_row_heatmap = "14_15"
flu_season_column_heatmap = "15_16"

rows_ha_mutations = []
rows_na_mutations = []
rows_tot_mutations = []

for sequence_name_row in key_flu_season[flu_season_row_heatmap]:
	row_ha_mutations = []
	row_na_mutations = []
	row_tot_mutations = []
	for sequence_name_column in key_flu_season[flu_season_column_heatmap]:
		seq1 = all_proteins_sequences[sequence_name_row]
		seq2 = all_proteins_sequences[sequence_name_column]
		
		ha_mutations = len(find_mutations(seq1["ha_sequence"], seq2["ha_sequence"]))
		na_mutations = len(find_mutations(seq1["na_sequence"], seq2["na_sequence"]))
		tot_mutations = ha_mutations + na_mutations
		
		row_ha_mutations.append(ha_mutations)
		row_na_mutations.append(na_mutations)
		row_tot_mutations.append(tot_mutations)
	rows_ha_mutations.append(row_ha_mutations)	
	rows_na_mutations.append(row_na_mutations)
	rows_tot_mutations.append(row_tot_mutations)

columns_name = [sequence_name_column for sequence_name_column in key_flu_season[flu_season_column_heatmap]]
rows_name = [sequence_name_row for sequence_name_row in key_flu_season[flu_season_row_heatmap]]
ha_mutations_df = pd.DataFrame(columns=columns_name)
na_mutations_df = pd.DataFrame(columns=columns_name)
tot_mutations_df = pd.DataFrame(columns=columns_name)

# Append each list in the list of lists as a new row to the dataframe
for row in rows_ha_mutations:
	ha_mutations_df = ha_mutations_df.append(pd.Series(row, index=columns_name), ignore_index=True)
for row in rows_na_mutations:
	na_mutations_df = na_mutations_df.append(pd.Series(row, index=columns_name), ignore_index=True)
for row in rows_tot_mutations:
	tot_mutations_df = tot_mutations_df.append(pd.Series(row, index=columns_name), ignore_index=True)

ha_mutations_df.index = rows_name
na_mutations_df.index = rows_name
tot_mutations_df.index = rows_name

ha_mutations_df = ha_mutations_df.apply(pd.to_numeric, errors='coerce')
na_mutations_df = na_mutations_df.apply(pd.to_numeric, errors='coerce')
tot_mutations_df = tot_mutations_df.apply(pd.to_numeric, errors='coerce')

# Create a heatmap using Seaborn
sns.clustermap(tot_mutations_df, cmap='coolwarm')
		