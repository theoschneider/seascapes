from libraries import filter_fasta

dico_fasta1 = {"seq1": "----ATGC-AT----",
               "seq2": "--ATATGC-AT---C"}

dico_fasta2 = {"seq1": "--AAATGC-AT----",
               "seq2": "--ATATGC-AT----"}

dico_fasta1, dico_fasta2 = filter_fasta(dico_fasta1, dico_fasta2)

print(dico_fasta1)
print(dico_fasta2)
