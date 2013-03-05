(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Basic Functions*)


(* ::Subsubtitle:: *)
(*Basic functions that provide usefull toosl for every other subpackage*)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Text:: *)
(*This sub package provides basic checks and commands to verfiy arguments of functions*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		13.11.2012*)
(*Last Changed: 	02.03.2013*)


BeginPackage["MPAWL`Basics`",{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by
Adriano Pascoletti, see
http://library.wolfram.com/infocenter/MathSource/7081/
*)
}];


(* ::Section:: *)
(*Global Function Declaration*)


isMatrixValid::usage="isMatrixValid[mM]

Check the validity of the matrix, i.e. whether mM is quadratic, regular and
integral.";


isMatrixValid::noMatrix = "The argument `1` is not a quadratic matrix";
isMatrixValid::noIntegerMatrix = "The argument `1` is not a an Integer matrix";
isMatrixValid::nonRegular = "The argument `1` is a nonregular matrix";


isDataValid::usage = "isDataValid[mM, data]

Check, whether the given data array is a set of points adressed in the cycles
of the (valid) matrix mM, i.e. it must be either of length m=Det[mM] or an 
array of Dimensions of the elementary divisors of mM that are greater than 1.

isDataValid[data,d, m, \[Epsilon]]
To avoid computation of the SNF, this method checks the same as above, where
d=Dimensions[mM][[1]], m=|Det[mM]| and \[Epsilon] denotes the elementary divisors of mM.";


isDataValid::ErrorInData = "The data is of dimension `1`, which is either not a set of `2` points in total or not of shape `3`.";


isIndexInRange::usage = "isIndexInRange[data,index]

Check whether the vector index can adress a value in data, i.e. whether each
entry of index corresponds to a valid entry in the corresponding dimension of
data. If the dimensions of data is greater than index, this adressing might
return a vector itself, even in that case this method checks dimensions.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)
Check \[Rule]  \!\(\*
StyleBox[\"\[OpenCurlyDoubleQuote]Positivity\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Negativity\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Nonzero\[CloseCurlyDoubleQuote]
	checks whether all indices are in Range and positive, negative or at least
	nonzero.";

isIndexInRange::tooManyIndices = "The index `1` has too many entries to adress a value or vector in the given data, which only consists of a `2`-dimensional Array.";
isIndexInRange::NoIndexVector = "The specified Index `1` is not a vector, which is neccessary to adress data in any array.";
isIndexInRange::NoValidCheck = "No valid Check Option was given.";

Options[isIndexInRange] = {Check -> "Positivity"};


saveCoefficients::usage = "saveCoefficients[waveletType, file, data]

Save the coefficients of data to the file. The first argument should
characterize the corresponding (scaling or wavelet) function completely, i.e.
contain the matrix mM and for the non Dirichlet case the dilation matrix mJ
and a characterization of g, to distinguish different shift invariant spaces.";


loadCoefficients::usage = "loadCoefficients[waveletType, file]

Load coefficients of data from file. The first argument should
characterize the corresponding (scaling or wavelet) function completely, i.e.
contain the matrix mM and for the non Dirichlet case the dilation matrix mJ
and a characterization of g, to check, whether it's the same shift invariant
space. If file does not contain that type of coefficients, \!\(\*StyleBox[\"Null\",FontWeight\[Rule]\"Bold\"]\) is returned.";


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


isMatrixValid[mM_] := Module[{d},
	d = Dimensions[mM][[1]];
	If [((Dimensions[Dimensions[mM]] != {2}) || (Dimensions[mM][[2]]!= d)),Message[isMatrixValid::noMatrix, MatrixForm[mM]];Return[False];];
	If[!ArrayQ[mM,_,IntegerQ],Message[isMatrixValid::noIntegerMatrix, MatrixForm[mM]];Return[False];];
	If[(Det[mM] == 0),Message[isMatrixValid::nonRegular,MatrixForm[mM]];Return[False];];
	Return[True];
];


isDataValid[mM_,data_] := Module[{d,m,epsilon,mE,dM},
	If[!isMatrixValid[mM], Return[$Failed]];
	d = Dimensions[mM][[1]];
	m = Det[mM];
	mE = IntegerSmithForm[mM,ExtendedForm-> False];
	dM = d+1- Position[Diagonal[mE], Min[Select[Diagonal[mE], # > 1 &]]][[1,1]];
	epsilon = Diagonal[mE][[d-dM+1;;d]];
	Return[isDataValid[d,m,epsilon,data]]
];


isDataValid[d_,m_,epsilon_,data_] := Module[{},
	If[((Dimensions[Dimensions[data]]=={1}) && (Dimensions[data][[1]]== m)),
		(* vectorial case *)Return[True]
	];
	If[(( Dimensions[Dimensions[data]] ==  Dimensions[epsilon]) && (Dimensions[data] == epsilon)),
		(*Matrix Case *)Return[True]
	];
(*No Valid Case*)
	Message[isDataValid::ErrorInData, Dimensions[data], m, epsilon];
	Return[False];
];


isIndexInRange[data_, index_, opts:OptionsPattern[Options[isIndexInRange]]] := Module[{indexnum,posInRange, positive, negInRange, negative},
If [Dimensions[Dimensions[data]][[1]] < Dimensions[index][[1]],
Message[isIndexInRange::tooManyIndices,MatrixForm[index, TableDirections-> Row],Dimensions[Dimensions[data]][[1]]];
Return[$Failed];
];
If[Dimensions[Dimensions[index]][[1]]!= 1,
Message[isIndexInRange::NoIndexVector,MatrixForm[index]];
Return[$Failed];
];
indexnum = Dimensions[index][[1]];
(* Check whether the first indexnum indices of data are big enough to handle
an adressing with the index *)
If [OptionValue[Check]=="Positivity",
Return[(( Abs[Dimensions[data][[1;;indexnum]]-index ]) == (Dimensions[data][[1;;indexnum]]-index ) 
(*all indices of index are <= dimensions[data] and all indices are positive*)
&&  (Abs[index-ConstantArray[1,indexnum]] == index-ConstantArray[1,indexnum] ))
];
];
If[OptionValue[Check]== "Negativity",
Return[(( Abs[Dimensions[data][[1;;indexnum]]+index ]) == (Dimensions[data][[1;;indexnum]]+index ) 
(*all indices of index are >=  dimensions[data] and all indices are negative *)
&&  (-Abs[index+ConstantArray[1,indexnum]] == index+ConstantArray[1,indexnum] ))
];
];
If[OptionValue[Check]== "Nonzero",
Return[(( Abs[Dimensions[data][[1;;indexnum]]-Abs[index]]) == (Dimensions[data][[1;;indexnum]]-Abs[index] ) 
(*all indices of index are in modulus <= dimensions[data] and nonzero *)
&&  (Abs[Abs[index]-ConstantArray[1,indexnum]] == Abs[index]-ConstantArray[1,indexnum] ))
];
];
Message[isIndexInRange::NoValidCheck];
Return[$Failed];
];


saveCoefficients[waveletType_,sFile_,data_] := Export[sFile,{waveletType,data},"List"];


loadCoefficients[waveletType_,sFile_] := Module[{temp},
	If[FileExistsQ[sFile],temp = Import[sFile,"List"],Return[$Failed]];
	If[(ToExpression[temp[[1]]]==waveletType),Return[ToExpression[temp[[2]]]]];
	Return[$Failed];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
