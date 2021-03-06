(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Functions providing functionality of translation invariant spaces*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		2012-11-13*)
(*Last Changed: 	2014-12-13*)


(* ::Subsubsection::Closed:: *)
(*License*)


(* ::Program:: *)
(*    This file is part of MPAWL.*)
(*  *)
(*      MPAWL is free software : you can redistribute it and/or modify*)
(*    it under the terms of the GNU General Public License as published by*)
(*    the Free Software Foundation, either version 3 of the License, or*)
(*    (at your option) any later version.*)
(*  *)
(*      MPAWL is distributed in the hope that it will be useful,*)
(*    but WITHOUT ANY WARRANTY; without even the implied warranty of*)
(*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the*)
(*    GNU General Public License for more details.*)
(*  *)
(*      You should have received a copy of the GNU General Public License*)
(*    along with the MPAWL. If not, see <http://www.gnu.org/licenses/>.*)


(* ::Subsection:: *)
(*Package Header*)


BeginPackage["MPAWL`TISpace`",
{
(*External dependencies*)
"SmithFormV6`", (* provided in this package, written by Adriano Pascoletti, see http://library.wolfram.com/infocenter/MathSource/7081/ *)
"MPAWL`Basics`",
"MPAWL`genSet`",
"MPAWL`Pattern`"
}];


(* ::Section:: *)
(*Global Function Declaration*)


(* ::Subsection:: *)
(*Bracket - Sums*)


computeBracketSums::usage = "computeBracketSums[data,originIndex,mM]

Compute the sum over the congruence classes h+mM^T*z, where h is from the
generating set and z runs through all integers addressing the values in data.

The result is given with respect of the coefficients of th generating set
basis, i.e. each h is decomposed with \!\(\*StyleBox[\"generatingSetBasisDecomp\", \"Code\"]\) and these
coefficients are used to address the sum of h in the result. Here, originIndex
denotes the index in data that corresponds to the origin.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM
	and the check, whether the Origin is in Range.

compute \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Bracket\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]absolute Squares\[CloseCurlyDoubleQuote]
	despite just summing up the entries, the second option sums up the
	absolute squares of the data entries.

Index \[Rule] \!\(\*StyleBox[\"None\",\nFontSlant\[Rule]\"Italic\"]\)
	If specified other than None, only this Index is computed and its value
	returned, provided it is in Range of the data.";


Options[computeBracketSums] = {MPAWL`Compute -> "Bracket", MPAWL`Validate -> True, Index -> None};


(* ::Subsection:: *)
(*getSpaceFromFourier & getFourierFromSpace*)


getSpaceFromFourier::usage = "getSpaceFromFourier[ckFun, ckSpace, originIndex, mM]

Let f be a function given by the Fourier coefficients ckFun and a function \[CurlyPhi]
also given by Fourier coefficients, ckSpace, both sharing the same originIndex.

Then this function returns the coefficients with respect o the matrix mM, such
that the Fourier transform of the result are the weights of the translates of \[CurlyPhi],
such that their sum yields f.

The array returned addresses each entry of the generating set with respect to
the generating set basis. If f is not in the translates of \[CurlyPhi], some entries
are either \[Infinity] or Immediate.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM
	and the check, whether the Origin is in Range of both Fourier coefficient
	arrays.";

Options[getSpaceFromFourier] = {MPAWL`Validate -> True};


getFourierFromSpace::usage ="getFourierFromSpace[coefficients, ckSpace, originIndex, mM]

The coefficients represent the Fourier transform of the weights which applied
to the translates --- with respect to mM - of a function \[CurlyPhi] result in a
function f and \[CurlyPhi] with its Fourier coefficients ckSpace (where originIndex
is the index representing the origin), this function reconstructs the Fourier
coefficients of f.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM
	and the check, whether the Origin is in Range of the Fourier coefficients.";


Options[getFourierFromSpace] = {MPAWL`Validate -> True};


getFourierFromSpace::OriginNotInRange = "The given Origin `2` is not in the Range of values given by data (Dimensions `1`)";
getFourierFromSpace::WrongHatAShape = "The coefficients are in wrong shape `1` with respect to the cycles of the generating set of the matrix `2`which are `3`.";


orthonormalizeTranslatesInSpace::usage="orthonormalizeTranslates[coefficients, mM, mJ]

Let V denote any mM-invariant space an orthonormal basis formed by the
translates of \[Xi]. The coefficients denote the Fourier transform
of the weights in the sum of translates of \[Xi] that represent a second
function \[CurlyPhi]\[Element]V. Suppose the translates of \[CurlyPhi] with respect to
mN = Inverse[mJ].mM are linear independent.

Then this function computes the Fourier transformed coefficients of the
weighted sum that represents \[CurlyPhi]\[Prime], which spans the same mN-invariant
space as \[CurlyPhi] but is orthonormalized.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM, mJ
	and mN.
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.";


orthonormalizeTranslatesInSpace::NotLinearIndependent = "The translates of the function represented by hata are not linear independent with respect to `1` as one of the divisors is zero";
orthonormalizeTranslatesInSpace::WrongDimensions = "The input given for `1` has The Dimensions `3`, but the matrix requires an input of dimension `2`.";


Options[orthonormalizeTranslatesInSpace] = {MPAWL`Debug -> "None", MPAWL`Validate -> True};



(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Helpers*)


CreateDirections::usage = "Generate all directions a sum can extend, i.e. all combinations of unit vectors and their negative values.";
CreateDirections[d_] := Module[{prevset},
  If[d==1,
    Return[{{1},{-1}}]
  ];
  prevset = CreateDirections[d-1];
  Return[Union[
    Append[#,-1] & /@prevset,
    Append[#,0] & /@prevset,
    Append[#,1] & /@prevset,
    {Append[ConstantArray[0,d-1],1]},
    {Append[ConstantArray[0,d-1],-1]}
    ]
  ];
]/;( IntegerQ[d] && d> 0);


(* ::Subsection:: *)
(*Bracket Sum functions*)


computeBracketSums[ckFun_,originIndex_,mM_,opts:OptionsPattern[]] := 
	localBracketSums[ckFun,originIndex,mM,OptionValue[Compute],
		OptionValue[MPAWL`Validate],OptionValue[Index]];


localBracketSums[ckFun_,originIndex_,mM_,cp_,True,index_] :=
	If[ (!isMatrixValid[mM]) || (!isIndexInRange[ckFun,originIndex]),
		Return[$Failed],
		Return[localBracketSums[ckFun,originIndex,mM,cp,False,index]]
	];


localBracketSums[ckFun_,originIndex_,mM_,cp_,False,None] :=
Module[{m,d,dM,epsilon,sums,sumsE,torigin,tmax,hM,dims,\[Epsilon],k,indicesM,indicesckF,mapIndices,rules},
    m = Det[mM];
    d = Dimensions[mM][[1]];
    epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-patternDimension[mM, validateMatrix -> False]+1;;d]];
    (* Result Indices *)
    If[cp == "Bracket",
        (* Map generates replacement rules, where the left hand size is an index from indicesM
            the right hand size are all indices that are mapped onto this index extracted from
                ckFun and summed up*)
    sumsE = ReplacePart[ConstantArray[0,epsilon],
			Normal[Total[Extract[ckFun,#]] &/@GroupBy[Tuples[Range[1,#] &/@ Dimensions[ckFun]],generatingSetBasisDecomp[#-originIndex,Transpose[mM]]+1 &]]]
    ];
    If[cp == "absolute Squares", (* same as above but all values are additionally Abs[*]^2*)
	    sumsE = ReplacePart[ConstantArray[0,epsilon],
			Normal[Total[Abs[Extract[ckFun,#]]^2] &/@GroupBy[Tuples[Range[1,#] &/@ Dimensions[ckFun]],generatingSetBasisDecomp[#-{2,2},Transpose[mM]]+1 &]]]
    ];
    Return[sumsE];
];


localBracketSums[ckFun_,originIndex_,mM_,cp_,False,index_?(VectorQ[#, IntegerQ] &)] :=
Module[{m,d,baseIndexDecomp,mapIndices,indicesckF,sumRes},
	m = Det[mM];
	d = Dimensions[mM][[1]];
	sumRes=0;
	baseIndexDecomp = generatingSetBasisDecomp[index-originIndex,Transpose[mM]];
	indicesckF = Tuples[Range[1,#] &/@ Dimensions[ckFun]];
	(* Mappings *)
	mapIndices = generatingSetBasisDecompVec[indicesckF-originIndex,Transpose[mM]]+1;
	If[cp == "Bracket",
		(* left hand size is an index from indicesM
			the right hand size are all indices that are mapped onto this index extracted from
				ckFun and summed up*)
		sumRes = Total[Extract[ckFun,Extract[indicesckF,Position[mapIndices,baseIndexDecomp]]]];
	];
	If[cp == "absolute Squares", (* same as above but all values are additionally Abs[*]^2*)
		sumRes = Total[Abs[Extract[ckFun,Extract[indicesckF,Position[mapIndices,baseIndexDecomp]]]]^2];
	];
	Return[sumRes];
]


(* All other v_ and index-Stuff *)
localBracketSums[ckFun_,originIndex_,mM_,cp_,v_,index_] := $Failed;


(* ::Subsection:: *)
(*getSpaceFromFourier & getFourierFromSpace*)


getSpaceFromFourier[ckFun_, ckSpace_, originIndex_, mM_, opts:OptionsPattern[]] :=
	localGetCoeffFromFun[ckFun, ckSpace, originIndex, mM, OptionValue[MPAWL`Validate]];


localGetCoeffFromFun[ckFun_, ckSpace_, originIndex_, mM_, True] :=
	If[ (!isMatrixValid[mM]) ||  (!isIndexInRange[ckSpace,originIndex])
		|| (!isIndexInRange[ckFun,originIndex]),
			Return[$Failed]
		,
			Return[localGetCoeffFromFun[ckFun,ckSpace,originIndex, mM, False]]
	];


localGetCoeffFromFun[ckFun_, ckSpace_, originIndex_, mM_, False] :=
Module[{hM,m,d, epsilon,tmax, torigin,checks,checksE,dims,actfactor,k,\[Epsilon]},
	hM = generatingSetBasis[Transpose[mM], Target -> "Symmetric", validateMatrix -> False];
	m = Det[mM];
	d = Dimensions[mM][[1]];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-patternDimension[mM, validateMatrix -> False]+1;;d]];
	tmax ={ Max[Ceiling[(Transpose[mM].#)[[1]] &/@{{-1/2,1/2},{1/2,1/2},{1/2,-1/2},{-1/2,-1/2}}]]+1, Max[Ceiling[(Transpose[mM].#)[[2]] &/@{{-1/2,1/2},{1/2,1/2},{1/2,-1/2},{-1/2,-1/2}}]]+1};
	torigin = tmax+1;
	checks = ConstantArray[Infinity,2tmax+1];
		(*Run through all elements of dataspace *)
	dims = Dimensions[ckSpace];
	Do[ 
		If[ckSpace[[Sequence@@(Table[\[Epsilon][k],{k,1,Length[dims]}])]] == 0,
			(*then the check-data must be zero and the value may than stay as it is *)
			If[ckFun[[Sequence@@(Table[\[Epsilon][k],{k,1,Length[dims]}])]] != 0,
				checks[[Sequence@@(modM[Table[\[Epsilon][k],{k,1,Length[dims]}]-originIndex, Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] = Indeterminate;
			];
		(* If Zero everything is okay and stays as it is (Infinity == arbitrary or the already given value *)
		, (*else isDataValid nonzero \[Rule] compute factor that would be needed for this entry *)
			If[!isIndexInRange[ckFun,Table[\[Epsilon][k],{k,1,Length[dims]}]],
				actfactor=0
			,
				actfactor = ckFun[[Sequence@@(Table[\[Epsilon][k],{k,1,Length[dims]}])]]/ckSpace[[Sequence@@(Table[\[Epsilon][k],{k,1,Length[dims]}])]];
			];
			(* No factor was set here before \[Rule] set*)
			If[checks[[Sequence@@(modM[Table[\[Epsilon][k],{k,1,Length[dims]}]-originIndex, Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] == Infinity,
				checks[[Sequence@@(modM[Table[\[Epsilon][k],{k,1,Length[dims]}]-originIndex, Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] = actfactor;
			];
			(* Different Factor was set \[Rule] Immediate *)
			If[checks[[Sequence@@(modM[Table[\[Epsilon][k],{k,1,Length[dims]}]-originIndex, Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] != actfactor,
				checks[[Sequence@@(modM[Table[\[Epsilon][k],{k,1,Length[dims]}]-originIndex, Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] = Indeterminate;
			];
		]; (*end dataspace nonzero*)
	,Evaluate[Sequence@@Table[{\[Epsilon][k],1,dims[[k]]},{k,1,Length[dims]}]]
	]; (*end do*)
	(*collect result in right cyrcles*)
	checksE = ConstantArray[0,epsilon];
	Do[
		checksE[[Sequence@@(Table[\[Epsilon][k],{k,1,Length[epsilon]}]+1)]]=checks[[ Sequence @@ (modM[Table[\[Epsilon][k],{k,1,Length[epsilon]}].hM,Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]];
		,Evaluate[Sequence@@Table[{\[Epsilon][k],0,epsilon[[k]]-1},{k,1,Length[epsilon]}]]
	];
	Return[checksE];
];


localGetCoeffFromFun[ckFun_, ckSpace_, originIndex_, mM_, v_] := $Failed;


getAllDirs[1] := {{1},{-1}};


getAllDirs[2] := {{1,1},{1,-1},{-1,1},{-1,-1}}


getAllDirs[d_] := Union[Table[Prepend[v,1],{v,#}],Table[Prepend[v,-1],{v,#}]] & /@{getAllDirs[d-1]};


getFourierFromSpace[coefficients_, ckSpace_, originIndex_, mM_, opts:OptionsPattern[]] :=
	localGetFunFromCoeff[coefficients, ckSpace, originIndex, mM, OptionValue[MPAWL`Validate]];


localGetFunFromCoeff[coefficients_, ckSpace_, originIndex_, mM_, True] :=
	If[!isMatrixValid[mM],
		$Failed,
		localGetFunFromCoeff[coefficients, ckSpace, originIndex, mM, False]
]	


localGetFunFromCoeff[coefficients_, ckSpace_, originIndex_, mM_, False] :=
	Module[{tempInd,tmax,torigin,coefficientsOI,d,epsilon,hM,result,dims,index,\[Epsilon],k},
	d = Dimensions[mM][[1]];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-patternDimension[mM, validateMatrix -> False]+1;;d]];
	hM = generatingSetBasis[Transpose[mM], Target -> "Symmetric", validateMatrix -> False];
	tmax = Table[Max[Ceiling[(1/2Transpose[mM].#)[[j]] &/@getAllDirs[d]]+1],{j,1,d}];
	torigin = tmax+1;
	coefficientsOI = ConstantArray[Infinity,2tmax+1];
	(*Pre rearrange*)
	Do[
		coefficientsOI[[Sequence @@ (modM[Table[\[Epsilon][k],{k,1,Length[epsilon]}].hM,Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]] = coefficients[[Sequence\[NonBreakingSpace]@@\[NonBreakingSpace](Table[\[Epsilon][k],{k,1,Length[epsilon]}]+1)]];
	,Evaluate[Sequence@@Table[{\[Epsilon][k],0,epsilon[[k]]-1},{k,1,Length[epsilon]}]]];
	dims = Dimensions[ckSpace];
	index = Table[\[Epsilon][k],{k,1,Length[dims]}];
    result = 
		Table[
			coefficientsOI[[ Sequence @@ (modM[index-originIndex,Transpose[mM],Target -> "Symmetric", validateMatrix -> False]+torigin)]]
				* ckSpace[[Sequence @@  (index)]]
		,Evaluate[Sequence@@Table[{\[Epsilon][k],1,dims[[k]]},{k,1,Length[dims]}]]];
	Return[result];
];


localGetFunFromCoeff[coefficients_, ckSpace_, originIndex_, mM_, v_] := $Failed;


orthonormalizeTranslatesInSpace[coeffs_,mM_,mJ_, opts:OptionsPattern[]] :=
	localOrthTInS[coeffs,mM,mJ,OptionValue[MPAWL`Debug],OptionValue[MPAWL`Validate]];


localOrthTInS[coeffs_,mM_,mJ_,db_,True] := Module[{d,dM,epsilon},
	If[!isMatrixValid[mM], Return[$Failed]];
	If[!isMatrixValid[mJ], Return[$Failed]];
	If[!isMatrixValid[Inverse[mJ].mM], Return[$Failed]];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]]; (*Cycle lengths*)
	If[Dimensions[coeffs] != epsilon,
		Message[orthonormalizeTranslatesInSpace::WrongDimensions,"the smaller space",MatrixForm[epsilon],MatrixForm[Dimensions[coeffs]]];
		Return[$Failed]
	];
	Return[localOrthTInS[coeffs,mM,mJ,db,True]];
];


localOrthTInS[coeffs_,mM_,mJ_,db_,False] := Module[{mN,NTg,InvNy,hN,d,dN,dM,\[Lambda]g,mu,mP,hataRes,epsilon,linIndep,t1,actBSq,k},
	d = Dimensions[mM][[1]];
	mN = Inverse[mJ].mM;
	dN = patternDimension[mN, validateMatrix -> False];
	dM = patternDimension[mM, validateMatrix -> False];
	If[StringCount[db,"Text"]>0,Print["Orthonormalizing coefficients..."]];
	NTg=Transpose[mN].Complement[generatingSet[Transpose[mJ]],{{0,0}}][[1]];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]]; (*Cycle lengths*)
	InvNy= Inverse[mN].((Complement[pattern[
							getPatternNormalform[mJ, validateMatrix -> False]
							,Target -> "Symmetric", validateMatrix -> False]
							,{{0,0}}])[[1]]);
	hN = generatingSetBasis[Transpose[mN], Target -> "Symmetric", validateMatrix -> False];
	\[Lambda]g = generatingSetBasisDecomp[NTg,Transpose[mM], Target -> "Symmetric", validateMatrix -> False];
	mu = Diagonal[IntegerSmithForm[mN, ExtendedForm-> False]][[d-dN+1;;d]]; (*Cycle lengths*)
	mP = Transpose[Table[
					generatingSetBasisDecomp[hN[[j]],Transpose[mM], Target -> "Symmetric", validateMatrix -> False]
					,{j,1,patternDimension[mN]}]];
	hataRes = ConstantArray[Infinity,epsilon];
	linIndep = True;
	t1 = AbsoluteTiming[
		Do[
			actBSq = Abs[coeffs[[Sequence@@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]]^2
					+ Abs[coeffs[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]] ]^2
					&[mP.Evaluate[Table[k[j],{j,1,dN}]]];
			If[actBSq == 0, Message[orthonormalizeTranslatesInSpace::NotLinearIndependent,MatrixForm[mN]];
			 				linIndep=False;Break[]
			];
			(hataRes[[Sequence@@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1) ]]= coeffs[[Sequence@@(modM[#,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]]*Sqrt[Abs[Det[mJ]]/actBSq])&[mP.Evaluate[Table[k[j],{j,1,dN}]]];
			(hataRes[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1) ]] = coeffs[[Sequence @@ (modM[#+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1) ]]*Sqrt[Abs[Det[mJ]]/actBSq])&[mP.Evaluate[Table[k[j],{j,1,dN}]]];
		,Evaluate[Sequence@@Table[{k[j],0,mu[[j]]-1},{j,1,dN}]]];][[1]];
	If[!linIndep,Return[$Failed]];
	If[StringCount[db,"Time"]>0,Print["Orthonormalization took ",t1," seconds."]];
	Return[hataRes];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
