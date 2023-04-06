from Bio import SeqIO
from Bio.Seq import Seq

list_of_positions_ha = [158, 213, 277]
list_of_positions_na = []

path_ha = "H3N2/selected_sequences_ha.fa"
path_na = "H3N2/selected_sequences_na.fa"

new_file = "H3N2/labels.fa"

labels_dict = {}
for seq_record in SeqIO.parse(path_ha, "fasta"):
	
	aminoacid_position_ha = "".join([seq_record.seq[position - 1] for position in list_of_positions_ha])
	labels_dict[seq_record.id] = aminoacid_position_ha

for seq_record in SeqIO.parse(path_na, "fasta"):
	
	aminoacid_position_na = "".join([seq_record.seq[position - 1] for position in list_of_positions_na])
	labels_dict[seq_record.id] += aminoacid_position_na

labels = []
for seq_record in SeqIO.parse(path_ha, "fasta"):
	seq_record.seq = Seq(labels_dict[seq_record.id])
	labels.append(seq_record)


print(labels)	
SeqIO.write(labels, new_file, "fasta")

