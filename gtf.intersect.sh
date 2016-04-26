gtf1=$1
gtf2=$2
name=$3

echo "Usage: sh gtf.intersect.sh compare.file compared.file"

bedtools intersect -a $1 -b $2 |awk '{for (i = 1; i <= NF; i++) {  if ($i ~ /gene_id|transcript_id|class_code/){  printf "%s %s", $i, $(i+1)  }  }  print ""}' |sort -u |sed 's/"//g' |sed 's/;/\t/g' |awk -v n=$name -F" " '{OFS="\t"; print $4,$2,$6,n}' > cuff.${name}.name
