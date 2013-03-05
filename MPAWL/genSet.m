(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*The generating set *)


(* ::Subsubtitle:: *)
(*Functions that provide access to generating sets, their creation and adressing*)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Text:: *)
(*This sub package provides basic checks and commands to verfiy arguments of functions*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		13.11.2012*)
(*Last Changed: 	02.03.2013*)


BeginPackage["MPAWL`genSet`",{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by
Adriano Pascoletti, see
http://library.wolfram.com/infocenter/MathSource/7081/
*),
"MPAWL`Basics`"
,
"MPAWL`Pattern`"
}
];


(* ::Section:: *)
(*Global Function Declaration*)


modM::usage = "modM[k, mM]

calculates the modulus of k with respect to the matrix mM, i.e. the vector h
in the unit cube such that k = h + mM*z.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Target \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Unit\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Symmetric\[CloseCurlyDoubleQuote]
	target domain of the modulus, eiter the unit cube or the unit cube shifted
	by -1/2.

validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM.";


Options[modM] := {Target \[Rule] "Unit", MPAWL`validateMatrix \[Rule] True};


(* ::Subsection:: *)
(*generating Set functions*)


generatingSetBasisDecomp::usage = "generatingSetBasisDecomp[k,mM]

For the standard Basis of the generating Set provided by
generatingSetBasis[mM] the (integer) Coefficients, that reconstruct x from the
basis (up to equivalence with respect to mod[k,mM].

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Target \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Unit\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Symmetric\[CloseCurlyDoubleQuote]
	target domain of the modulus, eiter the unit cube or the unit cube shifted
	by -1/2.

validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM.";


Options[generatingSetBasisDecomp] := {Target \[Rule] "Unit", MPAWL`validateMatrix \[Rule] True};


generatingSetBasis::usage = "generatingSetBasis[mM]

Returns patternDimension[mM] vectors, whose integral multiples (up to each
elementary divisor -1) reproduce the complete generating Set pattern. The
basis vectors are ordered with respect to nondecreasing cycle lenths
(elementary divisors).

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Target \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Unit\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Symmetric\[CloseCurlyDoubleQuote]
	target domain of the modulus, eiter the unit cube or the unit cube shifted
	by -1/2.

validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM.";


Options[generatingSetBasis] := {Target \[Rule] "Unit", MPAWL`validateMatrix \[Rule] True};



generatingSet::usage="generatingSet[mM]

Returns the set of integer vectors originalting from the corresponding
pattern[mM] by multiplying each element with mM.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Target \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Unit\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Symmetric\[CloseCurlyDoubleQuote]
	target domain of the modulus, eiter the unit cube or the unit cube shifted
	by -1/2.

validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM.";


Options[generatingSet] := {Target \[Rule] "Unit", MPAWL`validateMatrix \[Rule] True};



(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


modM[k_, mM_, opts:OptionsPattern[]] := 
  localModM[k, mM, OptionValue[Target], OptionValue[validateMatrix]];


(* usual Mod with MatrixCheck *)
localModM[k_, mM_, "Unit", True] := 
	If[isMatrixValid[mM], mM.(Mod[Inverse[mM].k, 1]), $Failed];


(* usual Mod without MatrixCheck *)
localModM[k_, mM_, "Unit", False] := mM.(Mod[Inverse[mM].k, 1]);


(* symmetric Mod with MatrixCheck *)
localModM[k_, mM_, "Symmetric", True] := 
	If[isMatrixValid[mM], mM.(Mod[Inverse[mM].k + 1/2, 1] - 1/2), $Failed];


(* symmetric ModM without MatrixCheck *)
localModM[k_, mM_, "Symmetric", False] :=  
	mM.(Mod[Inverse[mM].k + 1/2, 1] - 1/2)


(* All other versions with any third parameter lead to nothing, hence any other
OptionValue is not possibe*)
localModM[k_, mM_, s_, b_] := $Failed;


generatingSetBasisDecomp[k_, mM_, opts:OptionsPattern[]] := 
	localgenSetBDecomp[k,mM,OptionValue[Target],OptionValue[validateMatrix]];


localgenSetBDecomp[k_,mM_,t_,True] := 
	If[!isMatrixValid[mM],Return[$Failed],localgenSetBDecomp[k,mM,t,False]];


localgenSetBDecomp[k_,mM_,t_,False] := Module[{mE,mS,mP,d,dM,aBV},
	{mE,{mP,mS}} = IntegerSmithForm[Transpose[mM], ExtendedForm\[Rule] True];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM];
	aBV = Transpose[Inverse[mS]];
	Return[(modM[Inverse[aBV].k,mE, Target \[Rule] t])[[d-dM+1;;d]]];
];


localgenSetBDecomp[k_,mM_,t_,v_] := $Failed;


generatingSetBasis[mM_, opts:OptionsPattern[]] :=  
localGeneratingSetBasis[mM,OptionValue[Target],OptionValue[validateMatrix]];


localGeneratingSetBasis[mM_,t_,True] := 
	If[!isMatrixValid[mM],Return[$Failed],localGeneratingSetBasis[mM,t,False]];

localGeneratingSetBasis[mM_,t_,False] := Module[{mE,mP,mS,dM,d},
	{mE,{mP,mS}} = IntegerSmithForm[Transpose[mM], ExtendedForm \[Rule] True];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM];
	Return[Table[
		modM[Transpose[Inverse[mS]].UnitVector[d,d-dM+j],mM,
			Target \[Rule] t, validateMatrix \[Rule] False]
		, {j,1,dM}] ];
];


localgeneratingSetBasis[mM_,t_,v_] := $Failed;


generatingSet[mM_, opts:OptionsPattern[]] := 
	localGeneratingSet[mM, OptionValue[Target], OptionValue[validateMatrix]];


localGeneratingSet[mM_,t_,True] := If[!isMatrixValid[mM],
	Return[$Failed],localGeneratingSet[mM,t,False]
];


localGeneratingSet[mM_,t_,False] :=
(modM[mM.#,mM, Target \[Rule] t] & /@pattern[getPatternNormalform[mM], validateMatrix \[Rule] False]);


localGeneratingSet[mM_,t_,vM_] := $Failed;


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
