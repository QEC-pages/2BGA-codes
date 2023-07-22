#! This GAP program verifies parameteAOArs of the 2BGA codes constructed for the paper
#! Quantum two-block group algebra codes
#! by Hsiang-Ku Lin and Leonid P. Pryadko

## Using the script file, para.sh to find the parameter set. The 1st, 2nd, 3nd, 4th, 5th, and 6th columns are corresponding to the group size, group number, dimension, distance, and position sets of the group elements a and b excluding identity respectively.  
## example of a parameter set, 42   6   10   3   [ 3 ]   [ 4, 12, 24 ]
## Using the function, checkdata, to verifies the dimension and distance of the codes
## example, checkdata(SmallGroup(42,6),[3],[4,12,24]);
## Usinh the function, checkrank, to calculate the ranks of parity check matrices and the ranks of A, B, AB and delta_x and delta_z defined in the paper.
## example, checkrank(SmallGroup(42,6),[3],[4,12,24]);

# adjust the path in the second argument if needed
#SetPackagePath("QDistRnd","./QDistRnd" );
# package available from https://github.com/QEC-pages/QDistRnd/
LoadPackage("QDistRnd");

# adjust the path to run the gap program if needed 
#Read("MY_PATH/2BGA.gap");

# set up the code in binary field
F:=GF(2);

#LeftMat(elt,gg); Rightmat(elt,gg);
#Given a group, grp, and a sigle group element, elt, from GAP small group library, the function RightMat and LeftMat return the right matrix of elt, R(g) and L(g) respectively.
LeftMat:=function(elt,grp)
   local i,j,lis,n,miset,mi,mj,mip,veclist,vecperm,vec;
   lis:=Elements(grp);
   n:=Length(lis);
   miset:=[];
   veclist:=[];
   for j in [1..n] do
       mj:=lis[j];
       mi:=elt*mj;
       Append(miset,[mi]);
   od;
   for i in [1..n] do
       mip:=Position(miset,lis[i]);
       Add(veclist,mip);
   od;
   vecperm:=PermList(veclist);
   vec:=PermutationMat(vecperm,n,1);
   return vec;
end;



RightMat:=function(elt,grp)
    local i,j,lis,n,miset,mi,mj,mip,veclist,vecperm,vec;
    lis:=Elements(grp);
    n:=Length(lis);
    miset:=[];
    veclist:=[];
    for j in [1..n] do
        mj:=lis[j];
        mi:=mj*elt;
	Append(miset,[mi]);
    od;
    for i in [1..n] do
	mip:=Position(miset,lis[i]);
        Add(veclist,mip);
    od;
    vecperm:=PermList(veclist);
    vec:=PermutationMat(vecperm,n,1);
    return vec;
end;

# the functions Alcomb and Brcomb are used to calculate Eq. (36) in the paper.
# Alcomb([2,3,4],gg) Brcomb([5,6,15],gg)
# The 2nd argument, gg= SmallGroup(n1,n2). n1 and n2 are the values of the 1st and 2nd columns. The 1st arguemnt in the functions of Alcomb and Brcomb are the values of 5th and 6th columns from the parameters set.  
# Given a group, grp and a subset, elt1, the function Alcomb and Brcomb return the left and right matrices, A=I+L(g1)+L(g2)+... and B=I+R(g1)+R(g2)+... respectively.
# In gap, the whole group elements can be obtained from the command, Elements(grp). In the first argument of the function, [n1,n2,n3,..] represents the positions of those elements respectively. 

Alcomb:=function(elt1,grp)
  local elts,Al,id,n,list;
  elts:=Elements(grp);
  n:=Length(elts);
  id:=IdentityMat(n);
  Al:=(id+Sum(elt1,i->LeftMat(elts[i],grp)))*One(F);
  #Display(Al);
  return Al;
end;

Brcomb:=function(elt1,grp)
  local elts,Br,id,n,list;
  elts:=Elements(grp);
  n:=Length(elts);
  id:=IdentityMat(n);
  Br:=(id+Sum(elt1,i->RightMat(elts[i],grp)))*One(F);
  #Display(Br);
  return Br;
end;



# the function gxgzcomb return the parity matrices Hx and HzT for the Eq. (16) in the paper.
# gxgzcomb([2,3],[4,5],gg) or gxgzcomb([2],[4,5,15,18],gg)
# Given a group, grp and the subsets of elt1 [ga1,ga2,..] and elt2 [gb1,gb2,...], gxgzcomb retruns and parity matrices Hx and HzT.
# The 1st and 2nd arguments in the funtion are corresponding to the vaules of 5th and 6th columns in the parameter set (a and b defined in the paper). 

gxgzcomb:=function(elt1,elt2,grp)
   local list,elts,Al,Br,id,AA,BB,AAT,BBT,gx,gzT,n;
   elts:=Elements(grp);
   n:=Length(elts);
   id:=IdentityMat(n);
   #Print("id=",id,"\n");
   Al:=Sum(elt1,i->LeftMat(elts[i],grp));
   #Print("Al=",Al,"\n");
   Br:=Sum(elt2,i->RightMat(elts[i],grp));
   #Print("Br=",Br,"\n");
   AA:=TransposedMat(id+Al);
   #Print("AA=",AA,"\n");
   BB:=TransposedMat(id+Br);
   #Print("BB=",BB,"\n");
   BBT:=TransposedMat(id+TransposedMat(Br));
   AAT:=TransposedMat(id+TransposedMat(Al));
   gx:=TransposedMatMutable(Concatenation(AA,BB))*One(F);
   #Print("gx=","\n");
   #Display(gx);
   #Print("gx=",gx,"\n");
   gzT:=Concatenation(BBT,AAT)*One(F);
   #Print("gzT=","\n");
   #Display(gzT);
   list:=[gx,gzT];
   return list;
end;



# the function connected checks if a matrix is connected matrix or not.
# the function connected returns an empty list if a matrix is connected.

connected:=function(mat)
local i,len,check,checklist,che,eltset,eltset1,pos;
len:=Length(mat);
checklist:=[];
eltset:=[1..2*len];
check:=Positions(mat[1],Z(2)^0);
#Print("check=",check,"\n");
if not 1=check[1] then
SubtractSet(eltset,check);
Append(checklist,check);
#Print("checklist[1]=",checklist[1],"\n");
check:=Positions(mat[checklist[1]],Z(2)^0);
else
Append(checklist,check);
#Print("checklist=",checklist,"\n");
fi;
i:=1;

while Length(checklist)>0 and Length(eltset)>0 do
#Print("i=",i,"\n");
SubtractSet(eltset,check);
#Print("eltset=",eltset,"\n");
#Print("checklist=",checklist,"\n");
Remove(checklist,1);
if Length(checklist)>0 then
  che:=((checklist[1]-1) mod len)+1;
  #Print("first check list=",checklist[1],"; ","che=",che,"\n");
  check:=Positions(mat[che],Z(2)^0);
fi;
#Print("check=",check,"\n");
eltset1:=[1..2*len];
SubtractSet(eltset1,eltset);
SubtractSet(check,eltset1);
Append(checklist,check);
#Print("checklist=",checklist,"\n");
i:=i+1;
od;
#Print("checklist=",checklist,"\n");
#Print("eltset=",eltset,"  ",Length(eltset)=0,"\n");
return eltset;
end;



## the function checkrank check if rank(Hx)=rank(Hz)>rank(A)+rank(B)-rank(AB)
## checkrank(gg,[3],[3,5,8]), where gg:=SmallGroup(10,1),setA and setB are the subsets excluding Identity (a and b defined in the paper).
## the function returns set:=[len,ss,k,setA,setB,rankHx,rankHz,rankA,rankB,diff,diff]. The return list represents: the size of the group, len; the number of the group, ss; the dimension of the code, k; two elements a and b, setA and setB; the ranks of Hx, Hz, A, B, and AB, deltax and deltaz. 

checkrank:=function(grp,setA,setB)
local len,ss,gxz,k,rankHx,rankHz,al,br,rankA,rankB,rankAB,set,deltax,deltaz,diff;
len:=Size(grp);
ss:=IdSmallGroup(grp)[2];
Print("Abelian?",IsAbelian(grp),"\n");
gxz:=gxgzcomb(setA,setB,grp);
      if connected(gxz[1])=[] and connected(TransposedMat(gxz[2]))=[] then
        k:=2*len-RankMat(gxz[2])-RankMat(gxz[1]);
        if k>0 then
          rankHx:=RankMat(gxz[1]);
          rankHz:=RankMat(gxz[2]);
          al:=Alcomb(setA,grp)*One(F);
          br:=Brcomb(setB,grp)*One(F);
          rankA:=RankMat(al);
          rankB:=RankMat(br);
          rankAB:=RankMat(al*br);
	  diff:=rankA+rankB-rankAB;
	  deltax:=rankA+rankB-rankAB-rankHx;
	  deltaz:=rankA+rankB-rankAB-rankHz;
	            
          set:=[len,ss,k,setA,setB,rankHx,rankHz,rankA,rankB,rankAB,deltax,deltaz];
          #Print(set,"\n");
         
        fi;
      fi;
return set;

end;


# the function checkdata calculates the dimension and the distance of the codes given the group, grp, and group elements, a and b.
# the return list represts the group size id[1], the group number id[2], the dimension kk, the distance dd, the group eleemnts a and b and the last item is an empty list if the parity matrices Hx and HzT are both connected matrices. 

checkdata:=function(grp,elt1,elt2)
local gxz,rescon,id,len,kk,dd,final;
dd:=0;

id:=IdSmallGroup(grp);
len:=id[1];
gxz:=gxgzcomb(elt1,elt2,grp);
if connected(gxz[1])=[] and connected(TransposedMat(gxz[2]))=[] then rescon:="conntected"; else rescon:="not connected"; fi;
kk:=2*len-RankMat(gxz[2])-RankMat(gxz[1]);
if kk>0 then
        dd:=DistRandCSS(TransposedMat(gxz[2]),gxz[1],100000,1: field:=F,maxav:=10);
	fi;
final:=[id[1],id[2],kk,dd,elt1,elt2,rescon];
return final;
end;