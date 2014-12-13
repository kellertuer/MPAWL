(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Sampling*)


(* ::Subsubtitle:: *)
(*Functions that provide methods for sampling for functions and hence approximating by translates (via change of basis from interpolatory to (orthonormal) bases)*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		2012-11-14*)
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


BeginPackage["MPAWL`Sampling`",
{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by Adriano Pascoletti, see http://library.wolfram.com/infocenter/MathSource/7081/*),
"MPAWL`Basics`",
"MPAWL`Pattern`",
"MPAWL`genSet`",
"MPAWL`TISpace`",
"MPAWL`Transforms`"
}];


(* ::Section:: *)
(*Global Function Declaration*)


(* ::Text:: *)
(*Declaration of all functions that this package provides for global usage*)


sampleFunction::usage = "sampleFunction[mM, f]

Sample the function f on the pattern of mM on the 2\[Pi]-periodic torus.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Text&Time\[CloseCurlyDoubleQuote]
	activate text output, either just text, timings or both

File \[Rule] \!\(\*StyleBox[\"None\",\nFontSlant\[Rule]\"Italic\"]\)
	to specify a File where the resulting coefficients are stored and may be
	loaded the same way again. 
	
Method \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Point\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Point set\[CloseCurlyDoubleQuote] 
	How to perform the sampling, where usually the sampling is performed by calling
	f for each point individually. The other method creates a set X of all points
	and performs f[X], assuming, f can handle a set of points and returns a set
	of function values.
	
validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM.";


Options[sampleFunction] = {MPAWL`Debug -> "None", File -> None,
 	Method -> "Point", MPAWL`validateMatrix -> True};


changeBasis::usage = "changeBasis[mM,Coeffs, bracketSums]

Perform a change of Basis on the coefficients Coeffs, that are the DFT of some samples and
bracketSums represents the Bracket sums of a certain function, i.e. all are
nonzero such that the Lagrange function exists.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Input \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Frequency\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote]
	Domain of the discrete Input set. If \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] is given, a Fourier transform is
	performed before the change of basis

Output \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Frequency\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote]
	Domain of the discrete Output set. If \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] is given, a Fourier transform is
	performed after the change of basis

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
		whether to perform a check (via StyleBox[\"isMatrixValid\", \"Code\"]\)) on the matrix mM
		and the check, whether the Origin is in Range.

Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Text&Time\[CloseCurlyDoubleQuote]
	activate text output, either just text, timings or both.";


Options[changeBasis] = {MPAWL`Debug -> "None ", Input -> "Frequency",
 	Output -> "Frequency", MPAWL`Validate -> True};


changeBasis::WrongBasisShape = "The `1` are in wrong shape `2` with respect to the cycles of the generating set of the matrix `3` which are `4`.";


(* ::Section:: *)
(*Function Definition*)


(* ::Text:: *)
(*All function definitions are inside the private part, so only those functions are globally visible and thus usable, that are declared above via ::usage*)
(*And even their Error Messages are not visible to public*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Helpers*)


sampleFunction[mM_, f_, opts:OptionsPattern[]] := 
	localSample[mM,f,OptionValue[MPAWL`Debug],OptionValue[File],OptionValue[Method],
		OptionValue[validateMatrix]];


localSample[mM_,f_,db_,file_,method_,True] :=
	If[!isMatrixValid[mM],Return[$Failed], localSample[mM,f,db,file,method,False]];


localSample[mM_,f_,db_,None,"Point",False] := Module[{d,dM,epsilon,pMBasis,t,data,\[Epsilon]},
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	pMBasis = patternBasis[mM, validateMatrix -> False];
	If[StringCount[db,"Text"]>0,Print["Sampling the function..."]];
	t = AbsoluteTiming[
		data = Table[
			f[2\[Pi]*modM[Sum[\[Epsilon][j]*pMBasis[[j]],{j,1,dM}],IdentityMatrix[d],Target -> "Symmetric", validateMatrix -> False]]
			,Evaluate[Sequence@@Table[{\[Epsilon][j],0,epsilon[[j]]-1},{j,1,dM}]]
			];
		][[1]];
	If[StringCount[db,"Time"]>0,Print["Sampling took ",t," seconds."]];
	Return[data];
];


localSample[mM_,f_,db_,None,"Point set",False] := Module[{d,dM,epsilon,pMBasis,t,pointSet,values,\[Epsilon]},
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	pMBasis = patternBasis[mM, validateMatrix -> False];
	If[StringCount[db,"Text"]>0,Print["Sampling the function..."]];
	pointSet = Table[modM[
		Sum[\[Epsilon][j]*pMBasis[[j]], {j, 1, dM}], IdentityMatrix[d], 
	    Target -> "Symmetric", validateMatrix -> False]
					,Evaluate[Sequence@@
						Table[{\[Epsilon][j], 0, epsilon[[j]] - 1},{j,1,dM}]]
					];
	t = AbsoluteTiming[values = f[2 \[Pi] Flatten[pointSet, Range[dM]]];][[1]];
	If[StringCount[db,"Time"]>0,Print["Sampling took ",t," seconds."]];
	Return[Partition[values, epsilon[[2 ;; dM]]]];
];


localSample[mM_,f_,db_,file_String,method_,False] := Module[{data,t1},
	(* Try to load from File *)
	If[StringCount[db,"Text"]>0,Print["Loading data from file \[OpenCurlyDoubleQuote]",file,"\[CloseCurlyDoubleQuote]..."]];
	t1 = AbsoluteTiming[data = loadCoefficients[mM,file];][[1]];
	If[data==$Failed, (* loading did not work *)
		If[StringCount[db,"Text"]>0,Print["failed. Proceed with sampling..."]];
		data = localSample[mM,f,db,None,method,False];
		If[StringCount[db,"Text"]>0,Print["Saving data to file \[OpenCurlyDoubleQuote]",file,"\[CloseCurlyDoubleQuote]."]];
		saveCoefficients[mM,file,data];
	,
		If[StringCount[db,"Time"]>0,Print["done. Loading data took ",t1," seconds"]];	
	];
	Return[data];
];


localSample[mM_,f_,db_,file_,method_,v_] := $Failed; (* any other v that True/False and any other than String/None files *)


changeBasis[mM_,Coeffs_,BSums_, opts:OptionsPattern[]] := 
	localChangeBasis[mM,Coeffs,BSums,OptionValue[MPAWL`Debug], OptionValue[Input],
		OptionValue[Output],OptionValue[MPAWL`Validate]];


localChangeBasis[mM_,Coeffs_,BSums_,db_,in_,out_,True] := Module[{epsilon,d,dM},
	If[!isMatrixValid[mM], Return[$Failed]];
	d = Dimensions[mM][[1]];
	dM = patternDimension[mM, validateMatrix -> False];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]];
	If[Dimensions[BSums] != epsilon,
	 	Message[changeBasis::WrongBasisShape,"Bracket sums",MatrixForm[Dimensions[BSums]],MatrixForm[Transpose[mM]],MatrixForm[epsilon]];
	 	Return[$Failed];
	];
	If[Dimensions[Coeffs] != epsilon,
	 	Message[changeBasis::WrongBasisShape,"coefficients",MatrixForm[Dimensions[Coeffs]],MatrixForm[Transpose[mM]],MatrixForm[epsilon]];
	 	Return[$Failed];
	];
	Return[localChangeBasis[mM,Coeffs,BSums,db,in,out,False]];
];


localChangeBasis[mM_,Coeffs_,BSums_,db_,"Frequency","Frequency",False] := Module[{t1,c2},
	t1 = AbsoluteTiming[  c2 = 1/(Abs[Det[mM]]*BSums)*Coeffs;][[1]];
	If[StringCount[db,"Time"]>0,Print["The change of basis took ", t1, " seconds."]];
	Return[c2];
];


localChangeBasis[mM_,Coeffs_,BSums_,db_,"Time","Frequency",False] := Module[{t1,FTCoeffs},
	t1 = AbsoluteTiming[FTCoeffs = FourierTransformTorus[mM,Coeffs, MPAWL`Validate -> False];][[1]];
	If[StringCount[db,"Time"]>0,Print["The Fourier transform of the input took ", t1, " seconds."]];
	Return[localChangeBasis[mM,FTCoeffs,BSums,db,"Frequency","Frequency",False]]
];


localChangeBasis[mM_,Coeffs_,BSums_,db_,"Frequency","Time",False] := Module[{t1,FTCoeffs,data},
	FTCoeffs = localChangeBasis[mM,Coeffs,BSums,db,"Frequency","Frequency",False];
	t1 = AbsoluteTiming[data = FourierTransformTorus[mM,FTCoeffs, MPAWL`Validate -> False];][[1]];
	If[StringCount[db,"Time"]>0,Print["The Fourier transform of the result took ", t1, " seconds."]];
	Return[data];
];


localChangeBasis[mM_,Coeffs_,BSums_,db_,"Time","Time",False] := Module[{t1,FTCoeffs,data},
	t1 = AbsoluteTiming[FTCoeffs = FourierTransformTorus[mM,Coeffs, MPAWL`Validate -> False];][[1]];
	If[StringCount[db,"Time"]>0,Print["The Fourier transform of the input took ", t1, " seconds."]];
	FTCoeffs = localChangeBasis[mM,FTCoeffs,BSums,db,"Frequency","Frequency",False];
	t1 = AbsoluteTiming[data = FourierTransformTorus[mM,FTCoeffs, MPAWL`Validate -> False];][[1]];
	If[StringCount[db,"Time"]>0,Print["The Fourier transform of the result took ", t1, " seconds."]];
	Return[data];
];


localChangeBasis[mM_,Coeffs_,BSums_,db_,in_,out_,v_] := $Failed;


End[ ]


EndPackage[ ]
