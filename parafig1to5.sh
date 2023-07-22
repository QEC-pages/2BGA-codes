# The script file is used to collect data for Figs. 1 to 5 in the paper.                                          
# Enter the weight of left matrix, wa;                                                                            
# Enter the weight and dimension ranges of the codes in "for loop", wt and k.                                     
# The parameters are exported to the file "try.out".                                                              

name="nonabelian"
wa=2
tmp=tmp.dat
tmp1=tmp1.dat
tmp2=tmp2.dat

rm -f try.out

for k in 2
do
  for wt in {4,5,6,7,8}
  do
    rm -f $tmp
    rm -f $tmp1
    rm -f $tmp2
    rm -f try1.out
    unzip -p ${name}.zip diswt${name}_wt${wt}wtL${wa}_order*_k${k}.txt >> $tmp
    sort -k1,1n $tmp > $tmp1
    awk '($4 > dmin ) {dmin=$4 ; print}' $tmp1 >> $tmp2
    awk 'NR==1 {prev=$1; line=$0} $1!=prev {print line} {prev=$1; line=$0} END {if (NR>0) print}' $tmp2 > try1.out
    wb=$(expr ${wt} - ${wa})
    echo "## ${name} groups; Wa=${wa}; Wb=${wb}; k=${k}" >> try.out
    cat try1.out >> try.out
    echo >> try.out
    echo >> try.out
  done
done
