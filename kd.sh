# The script file is used to find parameters that satisfies the pattern
# kd > n, where the code length $n=2\ell$ is twice the group size $\ell$.
# Enter the type of group: "nonabelian" or "abelian";
# Enter the weight of the codes, 'wt'; the weight of left matrix, 'wa';
# Replacing 'diswt${name}_wt${wt}wtL${wa}_*.txt' with '*.txt' if interested in all codes.  
# The parameters are exported to the file "try.out".

name="nonabelian"
wt=8
wa=4

rm -f try.out

unzip -p ${name}.zip diswt${name}_wt${wt}wtL${wa}_*.txt | awk '$3 * $4 > 2 * $1 {print}' >> try.out


