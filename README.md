# 2BGA-codes

This repository contains the quantum error-correcting codes constructed by
Hsiang-Ku Lin and Leonid P. Pryadko for the paper [Quantum two-block group
algebra codes](http://arXiv.org/abs/2306.16400) (arXiv:2306.16400).

## Codes

The zip archives `nonabelian.zip` and `abelian.zip` contain the complete
information about all generated two-block group-algebra (2BGA) codes; please see
the paper for the details.  For example, the file
`diswtabelian_wt6wtL2_order50_k18.txt` from the file `abelian.zip` contains
information about connected 2BGA codes generated from abelian groups of order
`50`, with total row weight `6` and 1st group element of weight `2`, of
dimension `k=18`.  The file contains only one row
```
  50   5   18   2   [ 2 ]   [ 3, 4, 8 ]
```
where the parameters are:
- the group order $\ell=50$,
- group number `m=5` (in the `GAP` *Small Groups Library*),
- code dimension `k=18`,
- code distance `d=2`,
- and, finally, the lists of non-identity group element indices, in this case defining group algebra
  elements $a=1+g_2$ and $b=1+g_3+g_4+g_8$ (according to the default order of group elements as defined in `GAP`).

We used shell scripts to extract useful information from these `zip` archives.

## Scripts

- The shell script file `parafig1to5.sh` can be used to produce the data for Figs `1`
  to `5` in the paper.

- The script file `parafig6to9.sh` can be used to produce the data for Figs `6` to `9` in
  the paper.

- The script file `kd.sh` is used to find the codes whose parameters $[[n,k,d]]$ satisfy
  the inequality $kd > n$, where the code length $n=2\ell$ is twice the group size
  $\ell$. 

## `GAP` file
The file `2BGA.gap` contains two main `GAP` functions:  

* the function `checkdata()` to verify the distance and the dimension of
  the codes.  It takes three arguments: 
  - the group to use
  - the list of indices of non-identity group elements in $a$
  - the list of indices of non-identity group elements in $b$
     
* The function `checkrank()` to calculate the ranks of matrices `A`, `B`,
    `AB`, `Hx`, `Hz`, and the parameters $\delta_X$ and $\delta_Z$ defined in
    the paper.  It takes the same three arguments.

The indices of the group elements should correspond exactly to the default ordering of the group as it is defined in the `SmallGroups` library of `GAP`.
     
