(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Transforms*)


(* ::Subsubtitle:: *)
(*Functions that provide the multivariate Fourier and Wavelet Transform*)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Text:: *)
(*This sub package provides basic checks and commands to verfiy arguments of functions*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		13.11.2012*)
(*Last Changed: 	02.03.2013*)


BeginPackage["MPAWL`Transforms`",{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by
Adriano Pascoletti, see
http://library.wolfram.com/infocenter/MathSource/7081/
*),
"MPAWL`Basics`",
"MPAWL`Pattern`",
"MPAWL`genSet`"
}
];


(* ::Section:: *)
(*Global Function Declaration*)


(* ::Subsection:: *)
(*Introductionary functions*)


reshapeData::usage="reshapeData[M,data,direction]

Perform a reshape of data, where direction denotes
True: From vector to matrix
False: The other way around

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM
	and the check, whether the Origin is in Range.";


Options[reshapeData] := {MPAWL`Validate -> True};


(* ::Subsection:: *)
(*The Fourier Transform on the Torus*)


FourierTransformTorus::usage = "FourierTransformTorus[mM, b]

Perform the Fourier transform on the pattern with respect to mM. b is either a
vector of length m=|Det[mM]| or adressing the points with respect to the basis
of the pattern, i.e. the cycles having the length of the elementary divisors.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM
	and the check, whether the Origin is in Range.

Compute \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Numeric\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Exact\[CloseCurlyDoubleQuote]
	Providing numerical data, the Fourier method is used to perform the
	transform using FFT techniques. If all entries of mM and b are given exact,
	the \[OpenCurlyDoubleQuote]Exact\[CloseCurlyDoubleQuote] computation can be used to obtain the exact transform
";


Options[FourierTransformTorus] = {MPAWL`Validate -> True, MPAWL`Compute -> "Numeric"};


(* ::Subsection:: *)
(*Wavelet Transform*)


WaveletTransformTorus::usage="WaveletTransformTorus[mM, mJ, data, \[CurlyPhi]NCoeffs, \[Psi]NCoeffs]
WaveletTransformTorus[mM, mJ, ckData, ck\[CurlyPhi]M, ck\[CurlyPhi]N, ck\[Psi]N, originIndex]

For a decomposition of mM = mJ*mN, i.e. mN is integral and all three are regular,
and |Det[mJ]| = 2] this method performs a wavelet decomposition of the data. 

On the one hand, data and both subspaces are given as coefficients of the space
that is decomposed. Then the result consists of two arrays {dataS, dataW}
representing the two parts of data as coefficients of their corresponding spaces.

On the other hand, if data and all three spaces are given as Fourier
coefficients, the result are the two parts again as Fourier coefficients. This
second method is based on the first and hence slower.

\!\(\*StyleBox[\"Options\",\nFontWeight\[Rule]\"Bold\"]\)

Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Text&Time\[CloseCurlyDoubleQuote]
	activate text output, either just text, timings or both

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM
	and the check, whether the Origin is in Range.";


WaveletTransformTorus::WrongDimensions = "The input given for `1` has The Dimensions `3`, but the matrix requires an input of dimension `2`.";


Options[WaveletTransformTorus] = {MPAWL`MPAWL`Debug -> "None", MPAWL`Validate -> True};


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Introductionary functions*)


reshapeData[mM_,data_,direction_,opts:OptionsPattern[]]:=
	localReshapeData[mM,data,direction,OptionValue[MPAWL`Validate]];


localReshapeData[mM_, data_, direction_, True] :=
If[(!isMatrixValid[mM]) || (!isDataValid[mM,data]),
	Return[$Failed],
	Return[localReshapeData[mM, data, direction, False]]
];

localReshapeData[mM_, data_,direction_,False] := Module[{d,m,dM,epsilon,r},
	d = Dimensions[mM][[1]];
	m = Det[mM];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	If[((Dimensions[Dimensions[data]]=={1}) && (Dimensions[data][[1]]== m) && (epsilon[[dM]] != m)),
	If[!direction, Return[data]]; (*if we should reshape to vector and it is a vector*)
	Return[Table[
			data[[ 1+Sum[
			(*Sum over all adresses, where each cycleadress is spread over the
				following epsilons in product*)
			Product[epsilon[[k]],{k,j+1,Length[epsilon]}]*Evaluate[Subscript[r,j]],
			{j,1,Length[epsilon]}]
				]]
			,Evaluate[Sequence@@Table[{Subscript[r,i],0,epsilon[[i]]-1},{i,1,Length[epsilon]}]] (*run through all cycles*)
			]
		];
	, (* because of the validity of the data we are in the matrix case and only have to flatten if requested *)
	If[direction, Return[data], Return[Flatten[data]]];
	]; 
];


localReshapeData[mM_, data_, direction_, v_] := $Failed;


(* ::Subsection:: *)
(*Fourier Transform*)


FourierTransformTorus[mM_, b_,opts:OptionsPattern[]] :=
	localFTT[mM, b, OptionValue[MPAWL`Validate], OptionValue[Compute]];

localFTT[mM_, b_, True, cmp_] := Module[{},
	If[localReshapeData[mM, b, True, True] == $Failed, (*data or matrix not valid *)
		Return[$Failed]];
	If[ (cmp=="Numeric") 
		&& (!(And @@ (NumericQ[#] & /@Flatten[mM]))
			||!(And @@ (NumericQ[#] & /@Flatten[b]))
			), Return[$Failed]];
	If[ (cmp=="Exact") 
		&& (!(And @@ (NumericQ[#] & /@Flatten[b])) (* any number not exact? *)
			|| (Or @@ (InexactNumberQ[#] & /@Flatten[b]))  (* i.e. inexact number occuring?*)
			), Return[$Failed]];
	Return[localFTT[mM,b,False,cmp]] (*Call corresp. local Function w/out tests *)
];


localFTT[mM_, b_,False,"Exact"] :=
Module[{epsilon, m, d, dM,epsilonranges, internalb, resultb,flattenedinput},
	(*Abort if M is no fitting Matrix *)
	d = Dimensions[mM][[1]];
	m = Abs[Det[mM]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	flattenedinput = ((Dimensions[Dimensions[b]] == {1}) && (Dimensions[b][[1]] == m) && (epsilon[[dM]] != m));
	(*Flatten if neccessary and check dimensions else *)
	internalb = reshapeData[mM, b, True,MPAWL`Validate -> False];
	If[internalb == $Failed,Return[$Failed]];
	resultb = localFTTERec[epsilon, internalb];
	If[flattenedinput,resultb = Flatten[resultb]];
	Return[resultb]
];


localFTTERec[epsilon_,matrixb_] := Module[{e,copyb,temp,r,k},
If[Dimensions[epsilon] == {1},(*end of recursion *)
	e = Dimensions[matrixb][[1]];
	Return[Simplify[1/Sqrt[e]*Table[Sum[Exp[-2\[Pi] I j k / e]*matrixb[[j+1]], {j,0,e-1}], {k,0,e-1}]]]
];
(* All Epsilon-ranges on which we will perform an FFT, so excluding the last in epsilon *)
copyb = matrixb;
Do [
e = epsilon[[Length[epsilon]]];
copyb[[Sequence@@ Append[Table[Subscript[r,k]+1,{k,1,Length[epsilon]-1}],All]]] = 
Prepend[#[[ Length[#];;2;;-1]],First[#]] & [
	1/Sqrt[e]*Table[Sum[Exp[-2\[Pi] I j k / e]*copyb[[Sequence@@Append[Table[Subscript[r,k]+1,{k,1,Length[epsilon]-1}],j+1]]], {j,0,e-1}], {k,0,e-1}]
	];
,
Evaluate[Sequence@@ Table[{Subscript[r,i],0,epsilon[[i]]-1},{i,1,Length[epsilon]-1}]]
];
Do [
	copyb[[Sequence@@Append[Table[All,{k,1,Length[epsilon]-1}],e] ]] 
		= localFTTERec[epsilon[[1;;Length[epsilon]-1]], copyb[[Sequence@@Append[Table[All,{k,1,Length[epsilon]-1}],e] ]]];,
{e,1,epsilon[[Length[epsilon]]]}
];
Return[copyb];
];


localFTT[mM_, b_,False,"Numeric"] :=
Module[{epsilon, m, d, epsilonranges, internalb, resultb,flattenedinput,dM},
	d = Dimensions[mM][[1]];
	m = Abs[Det[mM]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	flattenedinput = ((Dimensions[Dimensions[b]] == {1}) && (Dimensions[b][[1]] == m) && (epsilon[[dM]] != m));
	(*Flatten if neccessary and check dimensions else *)
	internalb = reshapeData[mM, b, True, MPAWL`Validate -> False];
	resultb = localFTTNRec[epsilon, internalb];
	If[flattenedinput,resultb = Flatten[resultb]];
	Return[resultb]
];


localFTTNRec[epsilon_,matrixb_] := 
Module[{e,copyb,temp,r,k},
	If[Dimensions[epsilon] == {1},(*end of recursion *)
		Return[Fourier[matrixb]];
	];
	(* All Epsilon-ranges on which we will perform an FFT, so excluding the last in epsilon *)
	copyb = matrixb;
	Do [
		copyb[[Sequence@@ Append[Table[Subscript[r,k]+1,{k,1,Length[epsilon]-1}],All]]] = Fourier[copyb[[Sequence@@ Append[Table[Subscript[r,k]+1,{k,1,Length[epsilon]-1}],All]]]];
	,
	Evaluate[Sequence@@ Table[{Subscript[r,i],0,epsilon[[i]]-1},{i,1,Length[epsilon]-1}]]
	];
	Do [
		copyb[[Sequence@@Append[Table[All,{k,1,Length[epsilon]-1}],e] ]] 
		= localFTTNRec[epsilon[[1;;Length[epsilon]-1]], copyb[[Sequence@@Append[Table[All,{k,1,Length[epsilon]-1}],e] ]]];
	,{e,1,epsilon[[Length[epsilon]]]}
	];
	Return[copyb];
];


(* ::Subsection:: *)
(*Wavelet Transform*)


WaveletTransformTorus[mM_,mJ_,ckData_,ck\[CurlyPhi]M_,ck\[CurlyPhi]N_, ck\[Psi]N_,originIndex_,opts:OptionsPattern[]] :=
localWTT[mM,mJ,ckData,ck\[CurlyPhi]M,ck\[CurlyPhi]N, ck\[Psi]N,originIndex,OptionValue[MPAWL`Validate],OptionValue[MPAWL`Debug]] /;  (And @@ (ArrayQ[#,_,NumberQ] & /@ {ckData,ck\[CurlyPhi]M,ck\[CurlyPhi]N, ck\[Psi]N,originIndex}));


localWTT[mM_,mJ_,ckData_,ck\[CurlyPhi]M_,ck\[CurlyPhi]N_, ck\[Psi]N_,originIndex_,True,db_] :=
Module[{},
	If[!isMatrixValid[mM], Return[$Failed]];
	If[!isMatrixValid[mJ], Return[$Failed]];
	If[!isMatrixValid[Inverse[mJ].mM], Return[$Failed]];
	If[Abs[Det[mJ]]!= 2, Return[$Failed]];
	(* Check Origiin *)
	If[ (  (!isIndexInRange[ckData,originIndex])
		|| (!isIndexInRange[ck\[CurlyPhi]M,originIndex])
		|| (!isIndexInRange[ck\[CurlyPhi]N,originIndex])
		|| (!isIndexInRange[ck\[Psi]N,originIndex]) ),
		Return[$Failed]
	];
	Return[localWTT[mM,mJ,ckData,ck\[CurlyPhi]M,ck\[CurlyPhi]N,ck\[Psi]N,originIndex,False,db]];
];


localWTT[mM_,mJ_,ckData_,ck\[CurlyPhi]M_,ck\[CurlyPhi]N_, ck\[Psi]N_,originIndex_,True,db_] :=
Module[{sCoeffs,wCoeffs,data,dataS,dataW},
	sCoeffs = getSpaceFromFourier[ck\[CurlyPhi]N,ck\[CurlyPhi]M,originIndex,mM, MPAWL`Validate -> False];
	wCoeffs = getSpaceFromFourier[ck\[Psi]N,ck\[CurlyPhi]M,originIndex,mM, MPAWL`Validate -> False];
	data = getSpaceFromFourier[ckData,ck\[CurlyPhi]M,originIndex,mM, MPAWL`Validate -> False];
	{dataS,dataW}=localWTT[mM,mJ,data,sCoeffs,wCoeffs,False,db];
	Return[
		{getFourierFromSpace[dataS,ck\[CurlyPhi]N,originIndex,Inverse[mJ].mM, MPAWL`Validate -> False],
		getFourierFromSpace[dataW,ck\[Psi]N,originIndex,Inverse[mJ].mM, MPAWL`Validate -> False]}
	];
];


localWTT[mM_,mJ_,ckData_,ck\[CurlyPhi]M_,ck\[CurlyPhi]N_, ck\[Psi]N_,originIndex_,v_,db_] := $Failed /; (And @@ (ArrayQ[#,_,NumberQ] & /@ {ckData,ck\[CurlyPhi]M,ck\[CurlyPhi]N, ck\[Psi]N}));


WaveletTransformTorus[mM_,mJ_,data_,sCoeffs_,wCoeffs_,opts:OptionsPattern[]] :=
	localWTT[mM,mJ,data,sCoeffs,wCoeffs,OptionValue[MPAWL`Validate],OptionValue[MPAWL`Debug]];


localWTT[mM_,mJ_,data_,sCoeffs_,wCoeffs_,True,db_] := Module[{epsilon,d,dM},
	(*Ckecks*)
	If[!isMatrixValid[mM], Return[$Failed]];
	If[!isMatrixValid[Inverse[mJ].mM], Return[$Failed]];
	If[Abs[Det[mJ]]!= 2,Return[$Failed]];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	If[Dimensions[sCoeffs] != epsilon,
		Message[WaveletTransformTorus::WrongDimensions,"the smaller scaling space",MatrixForm[epsilon],MatrixForm[Dimensions[sCoeffs]]];
		Return[$Failed]
	];
	If[Dimensions[wCoeffs] != epsilon,
		Message[WaveletTransformTorus::WrongDimensions,"the wavelet space",MatrixForm[epsilon],MatrixForm[Dimensions[wCoeffs]]];
		Return[$Failed]
	];
	If[Dimensions[data] != epsilon,
		Message[WaveletTransformTorus::WrongDimensions, "the data",MatrixForm[epsilon],MatrixForm[Dimensions[data]]];
		Return[$Failed]
	];
	Return[localWTT[mM,mJ,data,sCoeffs,wCoeffs,False,db]];
];


(*FWT only working on Coefficients of the V_M space*)
localWTT[mM_,mJ_,data_,sCoeffs_,wCoeffs_,False,db_] := 
Module[{d,dM,dN,mN,hN,adN,adM,epsilon,mu,mP,NTg, dataS, dataW,\[Lambda]g,t1,k,j},
	(*Constants*)
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm -> False]][[d-dM+1;;d]];
	adM = Abs[Det[mM]];
	mN = Inverse[mJ].mM;
	dN = patternDimension[mN, validateMatrix -> False];
	mu = Diagonal[IntegerSmithForm[mN, ExtendedForm -> False]][[d-dN+1;;d]];
	hN = generatingSetBasis[Transpose[mN], Target -> "Symmetric", validateMatrix -> False];
	adN = Abs[Det[mN]];
	(* dyadic, hence NTg is unique *)
	NTg=Transpose[mN].(Complement[generatingSet[Transpose[mJ], validateMatrix -> False],{{0,0}}][[1]]);
	\[Lambda]g = generatingSetBasisDecomp[NTg,Transpose[mM], validateMatrix -> False];
	mP =Transpose[Table[generatingSetBasisDecomp[hN[[j]], Transpose[mM], validateMatrix -> False],{j,1,dN}]];
	(* 
			The Result 
	*)
	dataW = ConstantArray[0,mu];
	dataS = ConstantArray[0,mu];
	If[StringCount[db,"Text"]>0,Print["Performing the Wavelet transform..."]];
	t1 = AbsoluteTiming[
		Do[
			dataS[[Sequence@@( Table[Subscript[k, j],{j,1,dN}]+1)]] = 1/Sqrt[Abs[Det[mJ]]]*(
				Conjugate[sCoeffs[[Sequence@@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]]*data[[Sequence @@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]
					+ Conjugate[sCoeffs[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]]*data[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]&[mP.Evaluate[Table[Subscript[k, j],{j,1,dN}]]]);
			dataW[[Sequence@@ (Table[Subscript[k, j],{j,1,dN}]+1)]] = 1/Sqrt[Abs[Det[mJ]]]*(
				Conjugate[wCoeffs[[Sequence@@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]]*data[[Sequence @@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]
					+ Conjugate[wCoeffs[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]]*data[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]&[mP.Evaluate[Table[Subscript[k, j],{j,1,dN}]]]);
		,Evaluate[Sequence@@Table[{Subscript[k, j],0,mu[[j]]-1},{j,1,dN}]]];
	][[1]];
	If[StringCount[db,"Time"]>0,Print["Performing the Wavelet transform took ", t1, " seconds."]];
	Return[{dataS,dataW}];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
