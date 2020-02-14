# The script works on a file in association format:
# bash annotate.association_file.sh associations.txt > associations.annotated.txt
association_file=$1

annotation_file="annotation/SNP_annotations.hg19.txt"

header_file="header.annotated"

# Sort, join by snpID, sort by genomic position, output to stdout
cat $header_file <(join -t$'\t' <(sort -k1,1 $annotation_file ) <(sort -k1,1 $association_file ) | sort -k3,3V -k4,4n ) 
