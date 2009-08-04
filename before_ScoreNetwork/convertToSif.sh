awk '{ if (NF == 3) print $1 " pp " $2; else 
                    print $1 " pp " $3 }' $1 \
> $1.sif
