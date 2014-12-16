(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*The pattern*)


(* ::Subsubtitle:: *)
(*Functions that provide access to patterns, their creation and addressing*)


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


BeginPackage["MPAWL`Pattern`",{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by Adriano Pascoletti, seehttp://library.wolfram.com/infocenter/MathSource/7081/*),
"MPAWL`Basics`"
}
];


(* ::Section:: *)
(*Global Function Declaration*)


patternDimension::usage = "patternDimension[mM]

Returns the number of elementary divisors of mM, that are greater than 1.
This corresponds to the number of basis vectors of the pattern and hence the
dimension of the corresponding lattice.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM.";

Options[patternDimension] := {MPAWL`validateMatrix -> True};


(* ::Subsection:: *)
(*Pattern functions*)


patternBasis::usage = "patternBasis[mM]

Returns patternDimension[mM] vectors, whose integral multiples (up to each
elementary divisor -1) reproduce the complete pattern. They are ordered with
respect to nondecreasing cycle lengths (elementary divisors)

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)
validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM.";


Options[patternBasis] := {MPAWL`validateMatrix -> True};


pattern::usage = "pattern[mM]

generates the set of points inside the unit cube, whose multiplication with mM
results in an integral vector. The matrix mM must be in pattern normal form for
this fast algorithm to work, see \!\(\*
StyleBox[\"getPatternNromalform\", \"Code\"]\).

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)
Target \[Rule]  \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Unit\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Symmetric\[CloseCurlyDoubleQuote]
	target domain of the modulus, either the unit cube or the unit cube shifted
	by -1/2.
validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM.";


Options[pattern] := {Target -> "Unit", MPAWL`validateMatrix -> True};


getPatternNormalform::usage  = "getPatternNormalform[mM]

Return the normal form of the matrix mM, that generates the same pattern, i.e.
the corresponding matrix is in upper triangular form and the dominant value is
on the diagonal.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)
validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via \!\(\*StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM.";


Options[getPatternNormalform] := {MPAWL`validateMatrix -> True};


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


patternDimension[mM_, opts:OptionsPattern[]] := localPatternDimension[mM,OptionValue[validateMatrix]]


localPatternDimension[mM_,True] := If[!isMatrixValid[mM],Return[$Failed],Return[localPatternDimension[mM,False]]];


localPatternDimension[mM_,False] := Module[{mE,d,f},
	d = Dimensions[mM][[1]];
	mE = IntegerSmithForm[mM,ExtendedForm-> False];
	(* +1 due to indexing starting at 1 *)
	Return[d+1- Position[Diagonal[mE], Min[Select[Diagonal[mE], # > 1 &]]][[1,1]]];
];


(* ::Subsection:: *)
(*Pattern functions*)


patternBasis[mM_, opts : OptionsPattern[]] := localPatternBasis[mM,OptionValue[validateMatrix]]


localPatternBasis[mM_, True] := If[!isMatrixValid[mM],Return[$Failed],Return[localPatternBasis[mM,False]]];


localPatternBasis[mM_,False] := Module[{mE,mP,mS,d,dM,j},
	{mE,{mP,mS}} = IntegerSmithForm[mM, ExtendedForm-> True];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	Return[Table[mS.(1/(Diagonal[mE][[d-dM+j]])*UnitVector[d,d-dM+j]), {j,1,dM}] ]; 
]


pattern[mM_,opts:OptionsPattern[]] := localPattern[mM, OptionValue[Target], OptionValue[validateMatrix]];


localPattern[mM_, "Symmetric", True] := If[!isMatrixValid[mM],Return[$Failed],Return[localPattern[mM,"Symmetric",False]]];


localPattern[mM_,"Symmetric",False] := Module[{i,startset,min,max,step,d},
	d = Dimensions[mM][[1]];
	step  = 1/Abs [mM[[d,d]]];
	min = Ceiling[-(1/step)/2];
	max = Ceiling[(1/step)/2];
	startset = {};
	For[i=min, i< max, i++,
		startset = Union[startset,{i*step*UnitVector[d,d]}];
	];
	recursivePatternGeneration[mM,d-1,startset]
];


localPattern[mM_, "Unit", True] := If[!isMatrixValid[mM],Return[$Failed], Return[localPattern[mM,"Unit",False]]];


localPattern[mM_,"Unit",False] := (Mod[#,1] &/@localPattern[mM,"Symmetric",False]); (*yeah, a little lazy, might be not that efficient *)


localPattern[mM_,s_,t_] := $Failed;


(* little local helper to recursively construct the pattern *)
recursivePatternGeneration[mM_,actdim_,pElem_] := Module[{newset,d,min,max,step,tsum,i,y,x},
	(*Add the dimension actdim*)
	newset = {};
	d = Dimensions[mM][[1]];
	step = 1/Abs[mM[[actdim,actdim]]];
	Do[
		(*calc min and max including removal of shear*)
		tsum = Sum[mM[[actdim,i]]x[[i]],{i,actdim+1,d}];
		min = Ceiling[-1/(2 step)+tsum];
		max = Ceiling[1/(2 step)+tsum];
		For[i=min, i< max, i++,
			y=x+(i-tsum)*step *  UnitVector[d,actdim];
			newset = Union[newset,{y}];
		];
	,{x,pElem}];
	recursivePatternGeneration[mM,actdim-1,newset]
];
(*End of recursion - don't add any further elements, just return the elements, that were created up to now*)
recursivePatternGeneration[mM_,0,pElem_] := pElem


(* gcdOnRows[m_,ci,ri] does the euklidean algorithm with the elements a= m_[ci,ci] and b= m[ri,ci] and does the same manipulations also on the columns ci+1,...d of the matrix m. in the end
m_[ci,ci] is the gcd ofa and b} and the latter one is zero *)
gcdOnRows[m_,ci_,ri_] := 
	Module[{cm,fac},
		cm = m; (*computational matrix*)
		If[cm[[ri,ci]] != 0, (*only case where we have sth to do *)
			(*first assure that both a and b are nonnegative*)
			If[cm[[ci,ci]]<0,
				fac= If[cm[[ri,ci]]< 0,
						Ceiling[cm[[ci,ci]]/cm[[ri,ci]]], 
						Floor[cm[[ci,ci]]/cm[[ri,ci]]]
						];
				cm[[ci]] = cm[[ci]]-fac*cm[[ri]];
			];
			(*so cm[[ci,ci]] >= 0 now, if it is equal to zero, add (substract) one time row ri*)
			If[cm[[ci,ci] ]== 0,
				cm[[ci]] = cm[[ci]] + Sign[cm[[ri,ci]]]*cm[[ri]]
			];
			fac = Ceiling[cm[[ri,ci]]/cm[[ci,ci]]];
			If[fac==0,fac--];
			If[cm[[ri,ci]]<0,
				cm[[ri]] = cm[[ri]] - fac*cm[[ci]]
			];
			(* Finally start with the algorithm of euklid, slightly optimized slow version, because the fast one interchanges rows*)
			While[ cm[[ri,ci]]!= 0,
				If[Abs[cm[[ci,ci]]] > Abs[cm[[ri,ci]]],
					fac = Quotient[cm[[ci,ci]],cm[[ri,ci]]];
					If[Mod[cm[[ci,ci]],cm[[ri,ci]]]==0, 
						fac=fac-Sign[cm[[ri,ci]]]*Sign[cm[[ci,ci]]];
					]; (*we don't want to reach zero in this branch*)
					cm[[ci]] = cm[[ci]]-fac*cm[[ri]];
				,
					fac = Quotient[cm[[ri,ci]],cm[[ci,ci]]];
					cm[[ri]] = cm[[ri]] - fac*cm[[ci]];
				];
			];
		];
	(*return result*)
	cm
]


getPatternNormalform[mM_,opts:OptionsPattern[]] := localMat2PNF[mM,OptionValue[validateMatrix]];


localMat2PNF[mM_,True] :=Module[{},
	If[!isMatrixValid[mM],Return[$Failed]];
	Return[localMat2PNF[mM,False]];
]


localMat2PNF[mM_,False]:= Module[{d,col,row,row2,fac,cm},
	d=Dimensions[mM][[1]];
	cm = mM;
	For[col=1,col<d,col++,(*form each column into upper triangular*)
		For[row=col+1,row<=d,row++,(*handle each element below main diagonal*)
			cm = gcdOnRows[cm,col,row];
		];
	];
	(* now we have an upper triangular matrix, we just need a positiv diagonal*)
	For[row=1,row<= d, row++,
		If[cm[[row,row]]<0, cm[[row]] = -cm[[row]];]
	];
	(*compute the above values to be nonnegative in the matrix *)
	For[row=2,row <= d,row++,(*form each column into positiv values beupper triangular*)
		For[row2=row-1,row2>0,row2--,(*handle each element below main diagonal *)
			fac=0;
			(*if the entry is negative and smaller than row,row we stil need 1*)
			If[cm[[row2,row]]<0,fac = -Floor[cm[[row2,row]]/cm[[row,row]]]]; 
			If[cm[[row2,row]]>= cm[[row,row]],fac = -Floor[cm[[row2,row]]/cm[[row,row]]]];
			cm[[row2]] = cm[[row2]]+fac*cm[[row]];
		];
	];
	Return[cm];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
