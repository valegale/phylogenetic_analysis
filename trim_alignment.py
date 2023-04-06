from Bio import AlignIO

path_alignment = "H3N2/na_BD/msa_na_nucleotides.fa"
msa = list(AlignIO.parse("{}".format(path_alignment), "fasta"))[0]

	   
for col in range(msa.get_alignment_length()): #trim at ATG
	
	if msa[1, col] == "A" and msa[1, col + 1] == "T" and msa[1, col + 2] == "G": #searching in the first sequence
		position_start = col
		print("First full column is {}".format(col))
		break
	
for col in reversed(range(msa.get_alignment_length())): #trim at the last full column
		
   if not "-" in msa[:, col]:
	   position_end = col
	   print("Last full column is {}".format(col))
	   break
	

trimmed_msa = msa[:, position_start:position_end - 2] # trimming between ATG and removing the stop codon (position_end - 2)


# Write the trimmed alignment to a file
output_path = "H3N2/na_BD/trimmed_msa_na.fa"
with open(output_path, "w") as output_file:
    AlignIO.write(trimmed_msa, output_file, "fasta")
	