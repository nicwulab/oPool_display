mkdir cdhit
for x in 0.6 0.65 0.7 0.75 0.8 0.85
do
    cd-hit -i result/TableS1_filtered.fa -o cdhit/HA_Abs_paired_$x.fa -c $x -M 50000 -d 0 -T 60 -n 3.5 -aL 0.9 -s 0.95  -uS 0.2  -sc 1 -sf 1
done
