awk '{print $1 "\t1\t1\t"; print $2 "\t1\t1\t"}' $1 > $1.abundance
awk '{print $1 "\t1" }' $1 > $1.local
