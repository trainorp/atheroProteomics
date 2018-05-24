var="$(cat pepSeqsStr.txt)"

java -jar PeptideMatchCMD_1.0.jar -a query -i sprot_index -q $var -o out.txt 


