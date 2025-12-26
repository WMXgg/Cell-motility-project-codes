(* ::Package:: *)

lenGrid=40;(*size of the space*)
initialDensity=0.01;(*initial cell density of one cell type*)
initialBiomass=150;(*average initial biomass*)
biomassLimit=300;(*division biomass*)
initialBdis=UniformDistribution[{initialBiomass,biomassLimit}];
bionoise=30;(*maxmium noise for cell division*)
lenfun=2; (*the number of exchanging functions*)

move={1,1};(*move or not*)
movepro=0.9; (*moveprobability*)
speedM=1(*grid/per min; note that speedM*dt cannot exceed 1, otherwise the motility not be accurately captured*);
spdis=PoissonDistribution[speedM];
Ropro=1(*rotate probability controlled by the deviation of a Normaldistribution; A higher Ropro reflects a higher rotate probability*);
Rodis=NormalDistribution[0,Ropro];
NonRoPro=N@CDF[Rodis,\[Pi]/8]-N@CDF[Rodis,-(\[Pi]/8)](*probability of no rotation*);
qstres=100;
diffM={1,1,1,1,Sqrt[2],Sqrt[2],Sqrt[2],Sqrt[2]};(*weight for the diffusion\:987a\:5e8f\:4e3a\:5de6\:ff0c\:53f3\:ff0c\:4e0a\:ff0c\:4e0b\:ff0c\:5de6\:4e0a\:ff0c\:53f3\:4e0a\:ff0c\:5de6\:4e0b\:ff0c\:53f3\:4e0b*)

spreadCoefficient={1,1};(*diffusion coefficients of PG*)
rl={0.001,0.001}(*leakage rates*);
rl=(1/(lenGrid^2))*rl(*multipy the ratio between the intra and extra cellular volume *);
rt={100,100}(*take up rates*);
rt=(1/(lenGrid^2))*rt(*multipy the ratio between the intra and extra cellular volume *);
gr=0.1(*essential growth rate*);
kg={1,1};(*yield cofficient of amino acids for cell growth; this set was not used*)
Kg={1,1}(*half saturation constant for cell growth*);
Iin=2*Kg(*Initial intracellular concentrations of PG*);


T=10000;(*maxmium iteration number*)
record=10;(*frequency of data record*)
fill=0.95;(*stop once x% of grids are filled*)
dt=1(*min*); 


(*parameter sets generation*)
Paragenerator[]:=Module[{},
paraset={};
Do[
Di=SetPrecision[10^RandomReal[{-1,0}],4];(*diffusion coefficients of PG*)
rl=SetPrecision[10^RandomReal[{-5,-1}],4](*leakage rates*);
rl=(1/(lenGrid^2))*rl(*multipy the ratio between the intra and extra cellular volume *);
rt=SetPrecision[10^RandomReal[{-1,4}],4](*take up rates*);
rt=(1/(lenGrid^2))*rt(*multipy the ratio between the intra and extra cellular volume *);
gr=SetPrecision[10^RandomReal[{-3,0}],4](*essential growth rate*);
paraset=Join[paraset,{{Di,rl,rt,gr}}];
,{i,setnum}];
Export[FileNameJoin[{NotebookDirectory[],"paraset.csv"}],paraset];
];


(*initial distribution generation*)
Inidis[]:=Module[{},
dis={};
Do[
ranpos=RandomSample[DeleteDuplicates@RandomInteger[{1,lenGrid},{1000,2}],IntegerPart[lenfun*initialDensity*lenGrid^2]];
dis=Join[dis,{ranpos}];
,{i,disnum}];
Export[FileNameJoin[{NotebookDirectory[],"inidis-"<>ToString@initialDensity<>".csv"}],dis];
];

