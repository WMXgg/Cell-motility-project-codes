(* ::Package:: *)

(*============================
\:5b9a\:4e49\:ff1aINITIALIZE \:7528\:4e8e\:6a21\:578b\:521d\:59cb\:5316\:ff0c\:751f\:6210\:7f51\:683c\:77e9\:9635gridMatx\:548c\:7ec6\:80de\:77e9\:9635cellMatx
\:8f93\:5165\:ff1anull
\:8f93\:51fa\:ff1anull
\:66f4\:65b0\:ff1acellMatx1(cell type1), cellMatx1(cell type2),gridMatx
\:9884\:5148\:5b9a\:4e49\:ff1alenGrid\[LongDash]\[LongDash]\:7f51\:683c\:5927\:5c0f; initialDensity\[LongDash]\[LongDash]\:521d\:59cb\:7ec6\:80de\:5bc6\:5ea6; initialBilmass\[LongDash]\[LongDash]\:521d\:59cb\:751f\:7269\:91cf\:ff1blenFun\[LongDash]\[LongDash]\:529f\:80fd\:6570\:ff08\:5373\:516c\:5171\:7269\:8d28\:79cd\:7c7b\:ff09
(*CellMatx: 1:cell order, 2-3: position; 4: biomass; 5-6: genotype; 7-8: Intracelluar conc.; 9: motility or not?; 10: moving angle;11: growth rate*)
=============================*)
INITIALIZE[]:=Module[{temp,inipos},
num={};
cellData={};
gridData={};
countFunnumber={};
countFun={};
celltype=Select[Tuples[{0,1},lenfun],Total[#]==1&];

inipos=DeleteDuplicates[RandomInteger[{1,lenGrid},
{IntegerPart[lenfun*initialDensity*lenGrid^2],2}]];

len=(Length@inipos)/lenfun;
celltypelist=Flatten[Transpose[Table[celltype,IntegerPart[initialDensity*lenGrid^2]]],1];
movelist=Flatten[Transpose[Table[move,IntegerPart[initialDensity*lenGrid^2]]],1];

cellMatx=Table[Join[{i},inipos[[i]],{RandomVariate@initialBdis},celltypelist[[i]],Iin,{movelist[[i]]},{RandomVariate@Rodis},{0}],{i,Length@inipos}];

gridMatx=ConstantArray[Table[0,{j,lenfun}],{lenGrid,lenGrid}];
occupyfreq=Table[0,{i,1,lenGrid},{j,lenGrid}];
Celllin={};
];


(*============================
\:901a\:8fc7\:9605\:8bfb\:53c2\:6570\:96c6\:8868\:8fdb\:884c\:521d\:59cb\:5316
\:9884\:5148\:5b9a\:4e49\:ff1alenGrid\[LongDash]\[LongDash]\:7f51\:683c\:5927\:5c0f; initialDensity\[LongDash]\[LongDash]\:521d\:59cb\:7ec6\:80de\:5bc6\:5ea6; initialBilmass\[LongDash]\[LongDash]\:521d\:59cb\:751f\:7269\:91cf\:ff1blenFun\[LongDash]\[LongDash]\:529f\:80fd\:6570\:ff08\:5373\:516c\:5171\:7269\:8d28\:79cd\:7c7b\:ff09
(*CellMatx: 1:cell order, 2-3: position; 4: biomass; 5-6: genotype; 7-8: Intracelluar conc.; 9: motility or not?; 10: moving angle;11: growth rate*)
=============================*)
ReadInitialze[]:=Module[{temp,inipos},
num={};
cellData={};
gridData={};
countFunnumber={};
countFun={};
celltype=Select[Tuples[{0,1},lenfun],Total[#]==1&];

inipos=indis[[rep1]];

spreadCoefficient={parasets[[rep2,1]],parasets[[rep2,1]]};
rl={parasets[[rep2,2]],parasets[[rep2,2]]};
rt={parasets[[rep2,3]],parasets[[rep2,3]]};
gr=parasets[[rep2,4]];

len=(Length@inipos)/lenfun;
celltypelist=Flatten[Transpose[Table[celltype,IntegerPart[initialDensity*lenGrid^2]]],1];
movelist=Flatten[Transpose[Table[move,IntegerPart[initialDensity*lenGrid^2]]],1];

cellMatx=Table[Join[{i},inipos[[i]],{RandomVariate@initialBdis},celltypelist[[i]],Iin,{movelist[[i]]},{RandomVariate@Rodis},{0}],{i,Length@inipos}];

gridMatx=ConstantArray[Table[0,{j,lenfun}],{lenGrid,lenGrid}];
occupyfreq=Table[0,{i,1,lenGrid},{j,lenGrid}];
Celllin={};
];


(*movement function*)
Move[]:=Module[{rand,tcm,pos,ang,dis},
rand=Ordering[RandomReal[100,Length@cellMatx]];
Do[
tcm=cellMatx[[ce]];
pos=tcm[[2;;3]];
ang=tcm[[10]];

(*rotate?*)
Detang=N@RandomVariate@Rodis;
Detang=If[Detang<-Pi||Detang>Pi,Pi,Detang];
ang=ang+Detang;
(*move*)
ang=Mod[ang+Pi,2Pi]-Pi;
If[ang==-Pi,ang=ang+2*Pi];
dis=N@RandomVariate[spdis*dt];
If[RandomReal[]<movepro&&dis!=0,
movevector=FromPolarCoordinates[{dis,ang}];
newpos=Round[pos+movevector];
(*Are there any other cells on the way? If so, stop in the front of the blocker cell and reverse the movement.
If[IntersectingQ[Take[cellMatx,All,{2,3}],route]\[Equal]False,
route=Table[newpos+FromPolarCoordinates[{i,ang}],{i,1,dis}];
movevector=FromPolarCoordinates[{dis,ang}];
newpos=Round[pos+movevector];
newang==ang,
newpos=First[Nearest[Intersection[Take[cellMatx,All,{2,3}],route],pos]]-FromPolarCoordinates[{1,ang}]];*)
newang=ang,
newpos=pos;
newang=ang
];
newpos1=newpos;
newang1=newang;
(*Reverse when contact a boundary*)
If[newpos[[1]]==lenGrid||newpos[[2]]<lenGrid&&newpos[[2]]>1,newang=Pi+newang;newpos={newpos[[1]]-1,newpos[[2]]}];
If[newpos[[1]]==1||newpos[[2]]<lenGrid&&newpos[[2]]>1,newang=Pi+newang;newpos={newpos[[1]]+1,newpos[[2]]}];
If[newpos[[2]]==lenGrid||newpos[[1]]<lenGrid&&newpos[[1]]>1,newang=Pi+newang;newpos={newpos[[1]],newpos[[2]]-1}];
If[newpos[[2]]==1||newpos[[1]]<lenGrid&&newpos[[1]]>1,newang=Pi+newang;newpos={newpos[[1]],newpos[[2]]+1}];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
(*do not cross the boundary. The cell will be reverse when meet the boundary*)
If[newpos[[1]]<1,newpos={2-newpos[[1]],newpos[[2]]-2*(newpos[[1]]-1)*Tan[newang]};newang=Pi+newang,
If[newpos[[1]]>lenGrid,newpos={2*lenGrid-newpos[[1]],newpos[[2]]-2*(newpos[[1]]-lenGrid)*Tan[newang]};newang=Pi+newang,
If[newpos[[2]]<1,newpos={newpos[[1]]-2*(newpos[[2]]-1)*Cot[newang],2-newpos[[2]]};newang=Pi+newang,
If[newpos[[2]]>lenGrid,newpos={newpos[[1]]-2*(newpos[[2]]-lenGrid)*Cot[newang],2*lenGrid-newpos[[2]]};newang=Pi+newang
]]]];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
If[newpos[[1]]<1||newpos[[1]]>lenGrid||newpos[[2]]<1||newpos[[2]]>lenGrid,newpos=RandomInteger[{1,lenGrid},2]];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
(*Are there one cell already occupy the same position? If so, the cell will not move in this iteration*)
If[MemberQ[Table[cellMatx[[i,2;;3]],{i,Length@cellMatx}],Round@newpos]==True,
newpos=pos;
newang=ang
];

cellMatx[[ce,2;;3]]=newpos;
cellMatx[[ce,10]]=newang;

,{ce,rand}];
(*Are there several cells occupy the same position? If so, random delelte the overalp cells.*)
cellMatx=SortBy[DeleteDuplicatesBy[RandomSample[cellMatx,Length@cellMatx],#[[2;;3]]&],First];
];


(*movement function plus QS*)
MoveQS[]:=Module[{rand,tcm,pos,ang,dis,qs},
rand=Ordering[RandomReal[100,Length@cellMatx]];
cellpos=Take[cellMatx,All,{2,3}];
cellmatrix=Table[If[MemberQ[cellpos,{i,j}]==True,1,0],{i,1,lenGrid},{j,lenGrid}];
occupyfreq=occupyfreq+cellmatrix;
Do[
tcm=cellMatx[[ce]];
pos=tcm[[2;;3]];
ang=tcm[[10]];
qs=occupyfreq[[pos[[1]],pos[[2]]]];
(*if enough number of cells pass through the position, cell stop move*)
If[qs>qstres,Continue[],
(*rotate?*)
Detang=N@RandomVariate@Rodis;
Detang=If[Detang<-Pi||Detang>Pi,Pi,Detang];
ang=ang+Detang;
(*move*)
ang=Mod[ang+Pi,2Pi]-Pi;
If[ang==-Pi,ang=ang+2*Pi];
dis=N@RandomVariate[spdis*dt];
(*cells tends to stop moving in the center place*)
center={(lenGrid+1)/2,(lenGrid+1)/2};
centerdis=EuclideanDistance[pos,center];
movepro1=movepro*(centerdis/((lenGrid+1)/2));
(*movepro1=If[(pos[[1]]-((lenGrid+1)/2))^2+(pos[[2]]-((lenGrid+1)/2))^2<(lenGrid/4)^2,0,movepro];*)
If[RandomReal[]<movepro1&&dis!=0,
movevector=FromPolarCoordinates[{dis,ang}];
newpos=Round[pos+movevector];
(*Are there any other cells on the way? If so, stop in the front of the blocker cell and reverse the movement.
If[IntersectingQ[Take[cellMatx,All,{2,3}],route]\[Equal]False,
route=Table[newpos+FromPolarCoordinates[{i,ang}],{i,1,dis}];
movevector=FromPolarCoordinates[{dis,ang}];
newpos=Round[pos+movevector];
newang==ang,
newpos=First[Nearest[Intersection[Take[cellMatx,All,{2,3}],route],pos]]-FromPolarCoordinates[{1,ang}]];*)
newang=ang,
newpos=pos;
newang=ang
];
newpos1=newpos;
newang1=newang;
(*Reverse when contact a boundary*)
If[newpos[[1]]==lenGrid||newpos[[2]]<lenGrid&&newpos[[2]]>1,newang=Pi+newang;newpos={newpos[[1]]-1,newpos[[2]]}];
If[newpos[[1]]==1||newpos[[2]]<lenGrid&&newpos[[2]]>1,newang=Pi+newang;newpos={newpos[[1]]+1,newpos[[2]]}];
If[newpos[[2]]==lenGrid||newpos[[1]]<lenGrid&&newpos[[1]]>1,newang=Pi+newang;newpos={newpos[[1]],newpos[[2]]-1}];
If[newpos[[2]]==1||newpos[[1]]<lenGrid&&newpos[[1]]>1,newang=Pi+newang;newpos={newpos[[1]],newpos[[2]]+1}];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
(*do not cross the boundary. The cell will be reverse when meet the boundary*)
If[newpos[[1]]<1,newpos=N@{2-newpos[[1]],newpos[[2]]-2*(newpos[[1]]-1)*Tan[newang]};newang=Pi+newang,
If[newpos[[1]]>lenGrid,newpos=N@{2*lenGrid-newpos[[1]],newpos[[2]]-2*(newpos[[1]]-lenGrid)*Tan[newang]};newang=Pi+newang,
If[newpos[[2]]<1,newpos=N@{newpos[[1]]-2*(newpos[[2]]-1)*Cot[newang],2-newpos[[2]]};newang=Pi+newang,
If[newpos[[2]]>lenGrid,newpos=N@{newpos[[1]]-2*(newpos[[2]]-lenGrid)*Cot[newang],2*lenGrid-newpos[[2]]};newang=Pi+newang
]]]];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
If[newpos[[1]]<1||newpos[[1]]>lenGrid||newpos[[2]]<1||newpos[[2]]>lenGrid,newpos=RandomInteger[{1,lenGrid},2]];
newpos=Round@newpos;
newang=Mod[newang+Pi,2Pi]-Pi;
(*Are there one cell already occupy the same position? If so, the cell will not move in this iteration*)
If[MemberQ[Table[cellMatx[[i,2;;3]],{i,Length@cellMatx}],Round@newpos]==True,
newpos=pos;
newang=ang
];

cellMatx[[ce,2;;3]]=newpos;
cellMatx[[ce,10]]=newang;
];
,{ce,rand}];
(*Are there several cells occupy the same position? If so, random delelte the overalp cells.*)
cellMatx=SortBy[DeleteDuplicatesBy[RandomSample[cellMatx,Length@cellMatx],#[[2;;3]]&],First];
];


(*Production, Take up and Leakage*)
leakFun[]:=Module[{inner,out},
rand=Ordering[RandomReal[100,Length@cellMatx]];
Do[
tcm=cellMatx[[ce]];
pos=tcm[[2;;3]];
If[pos[[1]]<1||pos[[1]]>40||pos[[2]]<1||pos[[2]]>40,Continue[],
genotype=tcm[[5;;4+lenfun]];
biomass=tcm[[4]];
inner=tcm[[7;;6+lenfun]];
out=gridMatx[[pos[[1]],pos[[2]]]];
(*take up*)
takeup=Table[Min[rt*biomass*out[[i]]*dt,out[[i]]],{i,1,lenfun}];
inner=inner+takeup;
out=out-takeup;
(*leak*)
leak=Table[Min[rl*biomass*(inner[[i]]-out[[i]])*dt,inner[[i]]],{i,1,lenfun}];
gridMatx[[pos[[1]],pos[[2]]]]=out+leak;
afterdif=Table[Max[inner[[i]]-leak[[i]],0],{i,1,lenfun}];
(*producer maintain a constant concentration Iin of the amino acid that they 
can produce as a result of tight regulation of the amino acid production rates*)
afterproduce=Table[If[genotype[[i]]==1,Iin[[i]],afterdif[[i]]],{i,1,lenfun}];
cellMatx[[ce,7;;6+lenfun]]=afterproduce;
];
,{ce,rand}];
]


(*diffusion*)
ad=Drop[Tuples[{0,1,-1},2],1];
adjacentGrid[gri_]:=Module[{temp},
temp=Table[Join[gri+ad[[i]],{diffM[[i]]}],{i,Length@ad}];
temp=Select[temp,MemberQ[#,0]==False&&MemberQ[#,lenGrid+1]==False&];
temp
];
tempFun[in1_,in2_]:=Module[{t1,tp},
t1=adjacentGrid[{in1,in2}];
t1p=Take[t1,All,2];
t1M=Flatten@Take[t1,All,{3}];
tp=gridMatx[[in1,in2]];

tp=tp+spreadCoefficient*(t1M . Extract[gridMatx,t1p]-Total@diffM*tp)/Total@diffM;
tp
];
nutpubSpread[]:=Module[{},
rand=Ordering[RandomReal[100,lenGrid*lenGrid]];
posM=Flatten[Table[{i,j},{i,lenGrid},{j,lenGrid}],1];
gridMatx=Table[tempFun[posM[[rand[[i]],1]],posM[[rand[[i]],2]]],{i,lenGrid*lenGrid}];
gridMatx=Partition[gridMatx,lenGrid];
];


(*Cell growth*)
biopubChange[]:=Module[{rand,pos1,inner,tcm,B,DeltaB},
rand=Ordering[RandomReal[100,Length@cellMatx]];
Do[
tcm=cellMatx[[ce2]];
pos1=tcm[[2;;3]];
B=tcm[[4]];
genotype=tcm[[5;;4+lenfun]];
inner=tcm[[7;;6+lenfun]];
DeltaB=(Product[kg[[i]]*inner[[i]]/(Kg[[i]]+inner[[i]]),{i,Length@inner}])*gr*B*dt;
inner=Table[Max[0,N[inner[[i]]-(DeltaB/B)*inner[[i]]*(Table[1,lenfun]-genotype[[i]])]],{i,Length@inner}];(*amino acids are used*)
cellMatx[[ce2,4]]=B+DeltaB;
cellMatx[[ce2,11]]=DeltaB/(B*dt);
cellMatx[[ce2,7;;6+lenfun]]=inner;
,{ce2,rand}];
];


(*Cell division*)
cellReproduce[]:=Module[{rand,cellcopy,p1,rpos,del,copy},
rand=Ordering[RandomReal[100,Length@cellMatx]];
cellcopy=cellMatx;
numSplit=0;
celllin={};
Do[
If[cellMatx[[ce1,4]]<biomassLimit+RandomReal[{-bionoise,bionoise}],Continue[]];
numSplit++;
p1=cellcopy[[ce1,2;;3]];
adj=adjacentGrid[{p1[[1]],p1[[2]]}];
(*if the adjacentGrid is not full, random choose one free position to occupy*)
freeposlist=Select[adj,MemberQ[Table[cellcopy[[i,2;;3]],{i,Length@cellcopy}],#[[1;;2]]]==False&];
If[freeposlist!={},
rpos=RandomInteger[{1,Length@freeposlist}];
copy=cellcopy[[ce1]];
copy[[1]]=Max[Flatten@Take[cellcopy,All,{1}]]+1;
copy[[2;;3]]=freeposlist[[rpos]][[1;;2]];
copy[[4]]/=2;
Detang=RandomVariate@Rodis;
Detang=If[Detang<-Pi||Detang>Pi,Pi,Detang];
copy[[10]]=copy[[10]]+Detang;
copy[[10]]=Mod[copy[[10]]+Pi,2Pi]-Pi;
If[copy[[10]]==-Pi,copy[[10]]=copy[[10]]+2*Pi];(*random choose a moving ang*)
cellcopy[[ce1,4]]/=2;
cellcopy=Append[cellcopy,copy];
celllin=Join[celllin,{{1,cellcopy[[ce1,1]],copy[[1]]}}],
(*If[MemberQ[Table[cellcopy[[i,2;;3]],{i,Length@cellcopy}],adj[[r]][[1;;2]]]==False,
copy=cellMatx[[i]];
copy[[1]]=Length[cellcopy]+1;
copy[[2;;3]]=adj[[r]][[1;;2]];
copy[[4]]/=2;
cellcopy[[i,4]]/=2;
cellcopy=Append[cellcopy,copy];
Continue[]];*)
(*if the adjacentGrid is full, random choose one to compete for the space*)
If[RandomChoice[{0,1}]==1,
adj=adjacentGrid[{p1[[1]],p1[[2]]}];
rpos=RandomInteger[{1,Length@adj}];
droppos=Flatten@Position[Table[cellcopy[[i,2;;3]],{i,Length@cellcopy}],adj[[rpos]][[1;;2]]];
copy=cellcopy[[ce1]];
copy[[1]]=Max[Flatten@Take[cellcopy,All,{1}]]+1;
copy[[2;;3]]=adj[[rpos]][[1;;2]];
copy[[4]]/=2;
Detang=RandomVariate@Rodis;
Detang=If[Detang<-Pi||Detang>Pi,Pi,Detang];
copy[[10]]=copy[[10]]+Detang;
copy[[10]]=Mod[copy[[10]]+Pi,2Pi]-Pi;
If[copy[[10]]==-Pi,copy[[10]]=copy[[10]]+2*Pi];(*random choose a moving ang*)
cellcopy[[ce1,4]]/=2;
cellcopy=Append[cellcopy,copy];
celllin=Join[celllin,{{2,cellcopy[[ce1,1]],copy[[1]]}}];
cellcopy=Drop[cellcopy,droppos],
cellcopy[[ce1,4]]/=2;]
];
,{ce1,rand}];
cellMatx=cellcopy;

];


(*random order*)
randomorderQS[]:=Module[{},
ra=RandomSample[{1,2,3,4,5}];
Which[ra[[1]]==1,MoveQS[],
ra[[1]]==2,leakFun[],
ra[[1]]==3,nutpubSpread[],
ra[[1]]==4,biopubChange[],
ra[[1]]==5,cellReproduce[]];
Which[ra[[2]]==1,MoveQS[],
ra[[2]]==2,leakFun[],
ra[[2]]==3,nutpubSpread[],
ra[[2]]==4,biopubChange[],
ra[[2]]==5,cellReproduce[]];
Which[ra[[3]]==1,MoveQS[],
ra[[3]]==2,leakFun[],
ra[[3]]==3,nutpubSpread[],
ra[[3]]==4,biopubChange[],
ra[[3]]==5,cellReproduce[]];
Which[ra[[4]]==1,MoveQS[],
ra[[4]]==2,leakFun[],
ra[[4]]==3,nutpubSpread[],
ra[[4]]==4,biopubChange[],
ra[[4]]==5,cellReproduce[]];
Which[ra[[5]]==1,MoveQS[],
ra[[5]]==2,leakFun[],
ra[[5]]==3,nutpubSpread[],
ra[[5]]==4,biopubChange[],
ra[[5]]==5,cellReproduce[]];
];

randomordermove[]:=Module[{},
ra=RandomSample[{1,2,3,4,5}];
Which[ra[[1]]==1,Move[],
ra[[1]]==2,leakFun[],
ra[[1]]==3,nutpubSpread[],
ra[[1]]==4,biopubChange[],
ra[[1]]==5,cellReproduce[]];
Which[ra[[2]]==1,Move[],
ra[[2]]==2,leakFun[],
ra[[2]]==3,nutpubSpread[],
ra[[2]]==4,biopubChange[],
ra[[2]]==5,cellReproduce[]];
Which[ra[[3]]==1,Move[],
ra[[3]]==2,leakFun[],
ra[[3]]==3,nutpubSpread[],
ra[[3]]==4,biopubChange[],
ra[[3]]==5,cellReproduce[]];
Which[ra[[4]]==1,Move[],
ra[[4]]==2,leakFun[],
ra[[4]]==3,nutpubSpread[],
ra[[4]]==4,biopubChange[],
ra[[4]]==5,cellReproduce[]];
Which[ra[[5]]==1,Move[],
ra[[5]]==2,leakFun[],
ra[[5]]==3,nutpubSpread[],
ra[[5]]==4,biopubChange[],
ra[[5]]==5,cellReproduce[]];
];

randomordersessile[]:=Module[{},
ra=RandomSample[{2,3,4,5}];
Which[ra[[1]]==2,leakFun[],
ra[[1]]==3,nutpubSpread[],
ra[[1]]==4,biopubChange[],
ra[[1]]==5,cellReproduce[]];
Which[ra[[2]]==2,leakFun[],
ra[[2]]==3,nutpubSpread[],
ra[[2]]==4,biopubChange[],
ra[[2]]==5,cellReproduce[]];
Which[ra[[3]]==2,leakFun[],
ra[[3]]==3,nutpubSpread[],
ra[[3]]==4,biopubChange[],
ra[[3]]==5,cellReproduce[]];
Which[ra[[4]]==2,leakFun[],
ra[[4]]==3,nutpubSpread[],
ra[[4]]==4,biopubChange[],
ra[[4]]==5,cellReproduce[]];
];


(*Mainrunloop*)
MainloopQS[]:=Module[{},
INITIALIZE[];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
time1=AbsoluteTime[];
Do[randomorderQS[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"motile-QS","lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print[time2-time1]
];

Mainloopmotile[]:=Module[{},
INITIALIZE[];
Export[FileNameJoin[{NotebookDirectory[],"motile","celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile","griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
time1=AbsoluteTime[];
Do[randomordermove[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"motile","lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"motile","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"motile","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"motile","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print[time2-time1]
];

Mainloopsessile[]:=Module[{},
INITIALIZE[];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
time1=AbsoluteTime[];
Do[randomordersessile[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"nonmotile","lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print[time2-time1]
];


ParaloopQS[]:=Module[{},
indis=ToExpression@Import[FileNameJoin[{NotebookDirectory[],"inidis-"<>ToString@initialDensity<>".csv"}]];
parasets=Import[FileNameJoin[{NotebookDirectory[],"paraset.csv"}]];
time1=AbsoluteTime[];
Do[
Do[
ReadInitialze[];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Do[randomorderQS[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"motile-QS","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print["Distribution No. "<>ToString@rep1<>", Replicates No. "<>ToString@rep2<>" finished; Total time used: "<>ToString[N[(time2-time1)/60]]<>" min"];
,{rep2,paranum1,paranum2}];
,{rep1,disno1,disno2}];
];

Paraloopmotile[]:=Module[{},
indis=ToExpression@Import[FileNameJoin[{NotebookDirectory[],"inidis-"<>ToString@initialDensity<>".csv"}]];
parasets=Import[FileNameJoin[{NotebookDirectory[],"paraset.csv"}]];
time1=AbsoluteTime[];
Do[
Do[
ReadInitialze[];
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Do[randomordermove[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"motile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print["Distribution No. "<>ToString@rep1<>", Replicates No. "<>ToString@rep2<>" finished; Total time used: "<>ToString[N[(time2-time1)/60]]<>" min"];
,{rep2,paranum1,paranum2}];
,{rep1,disno1,disno2}];
];

Paraloopsessile[]:=Module[{},
indis=ToExpression@Import[FileNameJoin[{NotebookDirectory[],"inidis-"<>ToString@initialDensity<>".csv"}]];
parasets=Import[FileNameJoin[{NotebookDirectory[],"paraset.csv"}]];
time1=AbsoluteTime[];
Do[
Do[
ReadInitialze[];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[0]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[0]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Do[randomordersessile[];
If[celllin!={},Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"lindata","Lineage"<>ToString[t]<>".csv"}],celllin]];
If[Divisible[t,record]==True,
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
gridMatx1=Flatten[Flatten[Take[gridMatx,All,All,{1}],{3}],1];
gridMatx2=Flatten[Flatten[Take[gridMatx,All,All,{2}],{3}],1];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
];
If[Length@cellMatx>=fill*lenGrid*lenGrid,
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"celldata","CellMatrix"<>ToString[t]<>".csv"}],cellMatx];
Export[FileNameJoin[{NotebookDirectory[],"nonmotile","Dis"<>ToString@rep1,"data"<>ToString@rep2,"griddata","GridMatrix"<>ToString[t]<>".xlsx"}],"Sheets"->{"fun1"->gridMatx1,"fun2"->gridMatx2},"Rules"];
Break[]];
,{t,T}];
time2=AbsoluteTime[];
Print["Distribution No. "<>ToString@rep1<>", Replicates No. "<>ToString@rep2<>" finished; Total time used: "<>ToString[N[(time2-time1)/60]]<>" min"];
,{rep2,paranum1,paranum2}];
,{rep1,disno1,disno2}];
];
