java -jar PeptideMatchCMD_1.0.jar -a index -d uniprot_sprot_wvar.fasta -i sprot_index_wvar

var="$(cat pepSeqsStr.txt)"

java -jar PeptideMatchCMD_1.0.jar -a query -i sprot_index_wvar -q $var -o out.txt 


