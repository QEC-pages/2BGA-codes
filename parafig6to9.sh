# The script file is used to collect data for Fig. 6 and 7 in the paper.
# Enter the type of group: "nonabelian" or "abelian";
# Enter the weight of the codes, 'wt'; the weight of left matrix, 'wa';
# Enter the dimension range of the code, 'kmin' and 'kmax'.
# The parameters are exported to the file "try.out".

name="nonabelian"
wt=8
wa=2
wb=$(expr ${wt} - ${wa})
kmax=30
kmin=4
tmp=tmp.dat
tmp1=tmp1.dat
tmp2=tmp2.dat


rm -f try.out

echo "## ${name} groups; Wa=${wa}; Wb=${wb}" >> try.out

for ((k=${kmin}; k<=${kmax}; k=k+2)); do
  rm -f tmp.dat
  rm -f tmp1.dat
  rm -f tmp2.dat
  rm -f try1.out
  unzip -p ${name}.zip diswt${name}_wt${wt}wtL${wa}_order*_k${k}.txt >> $tmp
  sort -k1,1n $tmp > $tmp1
  awk '($4 > dmin ) {dmin=$4 ; print}' $tmp1 >> $tmp2
  awk 'NR==1 {prev=$1; line=$0} $1!=prev {print line} {prev=$1; line=$0} END {if       (NR>0) print}' $tmp2 >> try1.out
  echo "##k=${k}" >> try.out
  cat try1.out >> try.out
  echo >> try.out
  echo >> try.out
done

