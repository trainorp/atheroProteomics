file = '/home/patrick/gdrive/AthroProteomics/data/uniprot_sprot_wvar_digested_Mass400to6000.txt'
lines = []
with open(file) as f:
    for line in f:
        if ('HUMAN' in line) or ('P00722' in line):
            lines.append(line)

file2= '/home/patrick/gdrive/AthroProteomics/data/digest2.txt'
with open(file2, 'a') as f:
    for line in lines:
        f.write(line)