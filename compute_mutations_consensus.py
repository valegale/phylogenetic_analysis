from Bio import SeqIO

ha_path_consensus = "./../H3N2/removed_duplicates/ha_proteins/ha_consensus.fa"
na_path_consensus = "./../H3N2/removed_duplicates/na_proteins/na_consensus.fa"

first_bin = "06_07"
last_bin = "19_20"

sequences = {}
for seq_record in SeqIO.parse(ha_path_consensus, "fasta"):
	sequences[seq_record.id] = {}
	sequences[seq_record.id]["ha_sequence"] = seq_record.seq

for seq_record in SeqIO.parse(na_path_consensus, "fasta"):
	sequences[seq_record.id]["na_sequence"] = seq_record.seq

def find_mutations(seq1, seq2, protein = "ha"):
	return [protein + "_" + seq1[i - 1] + str(i) + seq2[i - 1] for i in range(1, len(seq1) + 1) if seq1[i - 1] != seq2[i - 1] and seq1[i - 1] != "-" and seq2[i - 1] != "-"]


def padding_zero(i):
	return str(i).zfill(2)
		
	
for i in range(int(first_bin[0:2]), int(last_bin[0:2])):
	bin1 = "{}_{}".format(padding_zero(i), padding_zero(i+1))
	bin2 = "{}_{}".format(padding_zero(i+1), padding_zero(i+2))
	seq1 = sequences[bin1]
	seq2 = sequences[bin2]
	
	ha_mutations = len(find_mutations(seq1["ha_sequence"], seq2["ha_sequence"]))
	na_mutations = len(find_mutations(seq1["na_sequence"], seq2["na_sequence"]))
	
	print ("From flu season {} to flu season {}".format(bin1, bin2))
	print ("HA: {} NA: {} Total: {}".format(ha_mutations, na_mutations, ha_mutations + na_mutations))
	
