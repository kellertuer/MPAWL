(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*de La Vall\[EAcute]e Poussin kernel as scaling functions and corresponding wavelets*)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		15.11.2012*)
(*Last Changed: 	04.03.2013*)


(* ::Subsubsection:: *)
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


BeginPackage["MPAWL`deLaValleePoussin`",
{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by
Adriano Pascoletti, see
http://library.wolfram.com/infocenter/MathSource/7081/
*),
"MPAWL`Basics`",
"MPAWL`Pattern`",
"MPAWL`genSet`",
"MPAWL`TISpace`",
"MPAWL`Visualization`",
"MPAWL`Transforms`"
}
];


(* ::Section:: *)
(*Global Function Declaration*)


dilationMatrix2D::usage = "dilationMatrix2D[letter]

represents the 2-dimensional dilation matrices available in the de la Vall\[EAcute]e
Poussin case. These are named by letters and an additional sign for some
matrices, in total these are: \[OpenCurlyDoubleQuote]X\[CloseCurlyDoubleQuote], \[OpenCurlyDoubleQuote]Y\[CloseCurlyDoubleQuote], \[OpenCurlyDoubleQuote]D\[CloseCurlyDoubleQuote], \[OpenCurlyDoubleQuote]T+\[CloseCurlyDoubleQuote], \[OpenCurlyDoubleQuote]T-\[CloseCurlyDoubleQuote], \[OpenCurlyDoubleQuote]S-\[CloseCurlyDoubleQuote] and \[OpenCurlyDoubleQuote]S+\[CloseCurlyDoubleQuote]";


pyramidFunction::usage = "pyramidFunction[\[Alpha],x]

The d-dimensional analog of the de la Vall\[EAcute]e Poussin mean
shrunken for generality on the symmetric unitcube, i.e. having
a support from -1/2-\[Alpha] to 1/2+\[Alpha] in each dimension.

\[Alpha] may be a nonnegative number less than 1/2 or an array of
d elements containing such numbers, where d is the length of x.";


(* ::Subsection:: *)
(*de La Vallee Poussin Kernels and subspaces*)


delaValleePoussinMean::usage = "delaValleePoussinMean[g,mM]

Generate the de la Vall\[EAcute]e Poussin Kernel \!\(\*SubsuperscriptBox[\(\[CurlyPhi]\), \(M\), \(\[EmptySet]\)]\) in Fourier coefficients
based on the function g and the matrix mM. While mM is an integral regular
matrix of dimension d*d, the function g has to fullfill the three properties:

1) nonnegativity
2) strictly positive on the shifted (symmetric) unit cube
3) For every x the sum over all integer shifts g(x+z) equals 1.

For simplicity g might be given as a nonnegative number t\[LessEqual]\!\(\*FractionBox[\(1\), \(2\)]\) which corresponds to
the pyramidal (tensor product of 1D de la Vall\[EAcute]e Poussin means) function, where the
support is -\!\(\*FractionBox[\(1\), \(2\)]\)-t to \!\(\*FractionBox[\(1\), \(2\)]\)+t in each dimension.
Or an array of such numbers, performing the same idea with different width in each direction,
hence the array must be the same dimension as each dimension of mM.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

BacketSums \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | True
	Compute the Bracket sums and return an array {ckV, BSV} consisting of the
	Fourier coefficients ckV and the Bracket Sums BSV of the kernel.
Orthonormalize \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	Perform an orthonormalization of the translates of the kernel with respect
	to mM.
validateMatrix \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM.
File \[Rule] \!\(\*StyleBox[\"None\",\nFontSlant\[Rule]\"Italic\"]\) | String | {String,String}
	save the Fourier coefficients of the kernel to a file. If the Bracket sums
	are computed too, there have to be two string, the first representing the
	kernel, the second the Bracket sums to be saved.
Support \[Rule] \!\(\*StyleBox[\"\!\(\*FractionBox[\(1\), \(2\)]\)\",\nFontSlant\[Rule]\"Italic\"]\) | p
	specifies the support of g, if it is given as a function, to extend the shifted unit cube
	by p in each dimension. This can also be given as a d-dimensional vector, where each entry
	is nonnegative and smaller than \!\(\*FractionBox[\(1\), \(2\)]\).
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.";


Options[delaValleePoussinMean] = {BracketSums -> False, File -> None,
 	MPAWL`Debug -> "None", Orthonormalize -> True, MPAWL`validateMatrix -> True, Support -> 0 };


delaValleePoussinMean::WrongFiles = "The specified File option `1` does not consist of two strings";
delaValleePoussinMean::WrongFile = "The specified File option `1` does not consist of one string";


delaValleePoussinSubspaces::usage="delaValleePoussinSubspaces[g,mM,mJ]

For any dyadic dilation matrix mJ and a function g fulfilling the properties:

1) nonnegativity
2) strictly positive on the shifted (symmetric) unit cube
3) For every x the sum over all integer shifts g(x+z) equals 1.

this method divides the mM-invariant space V (having dimension m=|Det[mM]) into
two orthogonal subspaces by returning two coefficient arrays, whose Fourier
transforms are weights to the sum of translates of any function, that spans V
with its translates. Both functions provide linear independence with respect to
mN = Inverse[mJ].mM

For simplicity g might be given as a nonnegative number t\[LessEqual]\!\(\*FractionBox[\(1\), \(2\)]\) which corresponds to
the pyramidal (tensor product of 1D de la Vall\[EAcute]e Poussin means) function, where the
support is -\!\(\*FractionBox[\(1\), \(2\)]\)-t to \!\(\*FractionBox[\(1\), \(2\)]\)+t in each dimension.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Orthonormalize \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	Perform an orthonormalization of the translates of both functions with
	respect to mN.
Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM, mJ
	and mN.
File \[Rule] \!\(\*StyleBox[\"None\",\nFontSlant\[Rule]\"Italic\"]\) | {String,String}
	save the coefficients of both functions to a file each.
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.";


delaValleePoussinSubspaces::WrongFiles = "The specified File option `1` does not consist of two strings";


Options[delaValleePoussinSubspaces] = {Orthonormalize -> True, 
	MPAWL`Debug -> "None", File -> None, MPAWL`Validate -> True};


(* ::Subsection:: *)
(*The Dirichlet Kernels*)


DirichletKernel::usage = "DirichletKernel[mM]

provides a dirichlet kernel, which is a special case of the de la Vall\[EAcute]e Poussin mean,
where g=0, hence pyramidalFunction[d,0] is used. Here, the same options as for
the de la Vallee Poussin mean apply.";


DirichletKernelSubspaces::usage = "DirichletKernelSubspaces[mM,mJ]

provides the dirichlet kernel based subspaces, i.e. this function is the special case
of the delaValleePoussinSubspaces, where g=0, hence pyramidalFunction[d,0] is used.
Here, the same options as for the de la Vallee Poussin Subspaces apply.";


(* ::Subsection:: *)
(*The Wavelet Decomposition*)


decomposeData2D::usage="decomposeData2D[g(Vec), JSet(Vec), mM, data]
decomposeData2D[l,g,JSet, mM, data, ck\[CurlyPhi]M]

Decompose the data given as coefficients with respect to ck\[CurlyPhi]\[CurlyPhi] this
method decomposes the data on different dialtionPaths given by JSet(Vec),
where the corresponding de la Vallee Poussin means (and wavelets) are given by
g(Vec) (see delaValleePoussinMean for restrictions on g). The dilation matrices
are given in terms of the characters used by dilationMatrix2D.

Both g and JSet might be given as vectors. If both are vectors, gVec has to be
one element longer than JSetVec. If one is given as a vector, the other is
expanded to a vector of corresponding length, having constant entries. This
length determines the number of levels of the wavelet decomposition, while g
also specifies additionally the starting space of mM.

If none of them is a vector, the second definition is needed to specify a level
depth of decomposition.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

ImagePrefix \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Img\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | String
	prefix that is used to save the images; the path in terms of the matrices
	is appended. This String may also contain Folders, that are relative to the
	Notebook.
ImageSuffix \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote].png\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | String
	suffix to determine the file format of the images.
DataPrefix \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]Data\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | String
	prefix that is used to save the data files; here an S or a W for the two
	subspaces is added and the decomposition path in terms of the matrices is
	appended. This String may also contain Folders, that are relative to the
	Notebook.
PlotResolution \[Rule] \!\(\*StyleBox[\"64\",\nFontSlant\[Rule]\"Italic\"]\) | nonnegative Integer
	resolution of the resulting Image generated by discreteFourierSeries
computeWavelet ->  \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to compute the wavelet part of the specific level and produce it's image
computeScale \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | True
	whether to compute the scale part of the specific level and produce it's image
Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
	whether to perform a check (via isMatrixValid[mM]) on the matrix mM and mJ(s)
	Validity of data. If some matrices mN in the decomposition are not valid,
	the algorithm continues on the other leaves.
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.

and any Option of discretePlotFourierSeries (including Plot and BarLegend), that
are passed on and any Option of Export (especially Resolution but not ImageSize)
Show (especially ImageSize). By the last two, Fonts are scaled, but a pixel size
of the image is not that easy to be computed.";


decomposeData2D::wrongGandmJs = "The specification of function(s) g(Vec) and the (set of) matrix-Char-Vectors mJSet(Vec) don't mach any valid combinaton";
decomposeData2D::wrongDimensionsM = "The matrix `1` is not a 2-dimensional matrix";
decomposeData2D::wrongDimensionsData = "A data array of dimension `1` is required, but the provided data is of dimension";


Options[decomposeData2D] = {MPAWL`Debug -> "None", ImagePrefix -> "Img-", ImageSuffix -> ".png", DataPrefix -> "Data",
PlotResolution -> 512, MPAWL`Validate -> True, computeWavelet -> True, computeScale -> False,
 (* passed to RasterizeFourierImage*) ReturnVal -> "ColorImage", ColorLegend -> True, Frame-> True,
 (* Standards for Export & Show*) ImageResolution -> 300, ImageSize -> 12*72/2.54, AllowRasterization -> Automatic};


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


dilationMatrix2D["X"] := {{2,0},{0,1}};
dilationMatrix2D["Y"] := {{1,0},{0,2}};
dilationMatrix2D["D"] := {{1,-1},{1,1}};
dilationMatrix2D["X+"] := {{2,0},{1,1}};
dilationMatrix2D["Y+"] := {{1,1},{0,2}};
dilationMatrix2D["X-"] := {{2,0},{-1,1}};
dilationMatrix2D["Y-"] := {{1,-1},{0,2}};


dilationMatrix2D["S+"] := {{2,0},{1,1}};
dilationMatrix2D["T+"] := {{1,1},{0,2}};
dilationMatrix2D["S-"] := {{2,0},{-1,1}};
dilationMatrix2D["T-"] := {{1,-1},{0,2}};


localdlVPMatChars = {"X","Y","D","S+","S-","T+","T-","X+","X-","Y+","Y-"};


pyramidFunction[\[Alpha]_,x_] := Which[ Abs[x] <= 1/2-\[Alpha],1, (1/2-\[Alpha]<Abs[x]) && (Abs[x]<1/2+\[Alpha]),(1/2+\[Alpha]-Abs[x])/(2\[Alpha]),True,0]/;((\[Alpha]>0)&&(\[Alpha]<= 1/2) && (NumberQ[x]))


pyramidFunction[0,x_] := Which[ Abs[x] < 1/2,1, Abs[x]== 1/2,1/2,True,0]/; (NumberQ[x])


pyramidFunction[\[Alpha]V_,xV_] := Product[pyramidFunction[\[Alpha]V[[j]],xV[[j]]],{j,1,Length[xV]}]/; (VectorQ[\[Alpha]V, (#>= 0) &&(#<= 1/2) &] && VectorQ[xV, NumberQ] && Length[\[Alpha]V]==Length[xV]) 


pyramidFunction[\[Alpha]_,xV_] := pyramidFunction[ConstantArray[\[Alpha],Length[xV]],xV] /;((\[Alpha]>= 0) &&(\[Alpha]<= 1/2) && VectorQ[xV, NumberQ]) 


delaValleePoussinMean[g_,mM_,opts:OptionsPattern[]] :=
	localdlVP[g,mM,
	OptionValue[MPAWL`Debug],OptionValue[Support],OptionValue[BracketSums],
	OptionValue[Orthonormalize],OptionValue[File],
	OptionValue[validateMatrix]];


getAllDirs[1] := {{1},{-1}};


getAllDirs[2] := {{1,1},{1,-1},{-1,1},{-1,-1}}


getAllDirs[d_] := Union[Table[Prepend[v,1],{v,#}],Table[Prepend[v,-1],{v,#}]] & /@{getAllDirs[d-1]};


(* Check pre *)
localdlVP[g_,mM_,db_,supp_,bs_,orth_,file_,True] :=
	If[!isMatrixValid[mM], $Failed,
		localdlVP[g,mM,db,supp,bs,orth,file,False]
	];


	(* If g is a number take the pyramidFunction and the old support *)
localdlVP[g_,mM_,db_,supp_,bs_,orth_,file_,False] :=
	localdlVP[pyramidFunction[ConstantArray[g,Dimensions[mM][[1]]],#] &,mM,db,ConstantArray[g,Dimensions[mM][[1]]],bs,orth,file,False] /; (NumberQ[g]);


	(* If g is an array of numbers, take pyramindFunction -> g also represents support *)
localdlVP[g_,mM_,db_,supp_,bs_,orth_,file_,False] :=
	localdlVP[pyramidFunction[g,#]&,mM,db,g,bs,orth,file,False] /; (ArrayQ[g,_,(#>=0) && (#<=1/2) &] && (Length[g]==Dimensions[mM][[1]]))


	(* Expand supp to a Vector *)
localdlVP[g_,mM_,db_,supp_,bs_,orth_,file_,False] :=
	localdlVP[g,mM,db,ConstantArray[supp,Dimensions[mM][[1]]],bs,orth,file,False] /; (NumberQ[supp]);


(* Just create localdlVP, no BS, no orth, no File *)
localdlVP[g_,mM_,db_,supp_,False,False,None,False] :=
Module[{d,adM,suppV, InvMt,t1,ck\[CurlyPhi]M,max,origin,x},
	d = Dimensions[mM][[1]];
	adM = Abs[Det[mM]];
	InvMt = Inverse[Transpose[mM]];
	max = Table[Max[Ceiling[((1/2+supp[[j]])Transpose[mM].#)[[j]] &/@getAllDirs[d]]]+1,{j,1,d}];
	origin = max+1;	
	If[StringCount[db,"Text"]>0,Print["Computing de la Vall\[EAcute]e Poussin coefficients..."]];
	t1 = Timing[ck\[CurlyPhi]M = Table[1/adM*g[InvMt.Table[Subscript[x,j],{j,1,d}]]
	,Evaluate[Sequence@@Table[{Subscript[x, j],-max[[j]],max[[j]]},{j,1,d}]]];][[1]];
	If[StringCount[db,"Time"]>0,Print["Creating the de La Vall\[EAcute]e Poussin kernel took ",t1," seconds."]];		
	Return[ck\[CurlyPhi]M];
]/; (VectorQ[supp,NumberQ] && Length[supp] == Dimensions[mM][[1]]);


(* Orthonormalize case *)
localdlVP[g_,mM_,db_,supp_,False,True,None,False] :=
Module[{\[CurlyPhi]MBSq,t,origin,max,ck\[CurlyPhi]M,adM,d},
	ck\[CurlyPhi]M = localdlVP[g,mM,db,supp,False,False,None,False];
	d = Dimensions[mM][[1]];
	max = Table[Max[Ceiling[((1+2supp[[j]])Transpose[mM].#)[[j]] &/@{{-1/2,1/2},{1/2,1/2},{1/2,-1/2},{-1/2,-1/2}}]]+1,{j,1,d}];
	origin = max+1;
	adM = Abs[Det[mM]];
	If[StringCount[db,"Text"]>0,Print["Orthonormalization..."]];
	t = Timing[
		\[CurlyPhi]MBSq = N[computeBracketSums[ck\[CurlyPhi]M,origin,mM,Compute -> "absolute Squares", MPAWL`Validate -> False]];
		ck\[CurlyPhi]M = getFourierFromSpace[1/(Sqrt[adM*\[CurlyPhi]MBSq]),ck\[CurlyPhi]M,origin,mM, MPAWL`Validate-> False]
	][[1]];
	If[StringCount[db,"Time"]>0,Print["Orthonormalization took ",t," seconds."]];		
	Return[ck\[CurlyPhi]M];
];


(* With Bracket Sums (with or without Orth, don't care)*)
localdlVP[g_,mM_,db_,supp_,True,orth_,None,False] :=
Module[{ck\[CurlyPhi]M, \[CurlyPhi]BS,max,origin,d,t1},
	d = Dimensions[mM][[1]];
	ck\[CurlyPhi]M = localdlVP[g,mM,db,supp,False,orth,None,False];
	max = Table[Max[Ceiling[((1+2supp[[j]])Transpose[mM].#)[[j]] &/@{{-1/2,1/2},{1/2,1/2},{1/2,-1/2},{-1/2,-1/2}}]]+1,{j,1,d}];
	origin = max+1;	
	If[StringCount[db,"Text"]>0,Print["Computing Bracket sums to obtain the basis transform coefficients..."]];
	{t1,\[CurlyPhi]BS} = AbsoluteTiming[computeBracketSums[ck\[CurlyPhi]M,origin,mM,Compute -> "Bracket", MPAWL`Validate -> False]];
	If[StringCount[db,"Time"]>0,Print["Computing the Bracket sums took ", t1, " seconds."]];
	Return[{ck\[CurlyPhi]M,\[CurlyPhi]BS}];
];


(* File Cases *)
(* a) One File only occurs without BS *)
localdlVP[g_,mM_,db_,supp_,False,orth_,file_String,False] :=
Module[{ck\[CurlyPhi]M,t1},
	If[StringCount[db,"Text"]>0,Print["Loading de la Vall\[EAcute]e Poussin mean from file \[OpenCurlyDoubleQuote]",file,"\[CloseCurlyDoubleQuote]..."]];
	{t1,ck\[CurlyPhi]M} = AbsoluteTiming[loadCoefficients[Join[mM,supp],file]];
	If[ck\[CurlyPhi]M == $Failed,
		If[StringCount[db,"Text"]>0,Print["failed. "]];
		ck\[CurlyPhi]M = localdlVP[g,mM,db,supp,False,orth,None,False];
		saveCoefficients[Join[mM,supp],file,ck\[CurlyPhi]M];
	, (* Successful: *)
	If[StringCount[db,"Time"]>0,Print["Loading coefficients took ",t1," seconds."]];	
	];
	Return[ck\[CurlyPhi]M];
];


(* BS but not one files or None *)
localdlVP[g_,mM_,db_,supp_,False,orth_,s_,False] := Module[{},
	Message[delaValleePoussinMean::WrongFile,s]; Return[$Failed];
];


(* b) Two Files only occur with BS *)
localdlVP[g_,mM_,db_,supp_,True,orth_,{ckFile_String,bsFile_String},False] :=
Module[{ck\[CurlyPhi]M,\[CurlyPhi]MBS,t1,t2,max,origin,d},
	(* Obtain ck\[CurlyPhi] via previous case *)
	ck\[CurlyPhi]M = localdlVP[g,mM,db,supp,False,orth,ckFile,False];
	If[StringCount[db,"Text"]>0,Print["Loading Bracket sums from file \[OpenCurlyDoubleQuote]",bsFile,"\[CloseCurlyDoubleQuote]..."]];
	{t1,\[CurlyPhi]MBS} = AbsoluteTiming[loadCoefficients[Join[mM,supp],bsFile]];
	If[\[CurlyPhi]MBS == $Failed,
		If[StringCount[db,"Text"]>0,Print["failed. Computing Bracket sums to obtain the basis transform coefficients..."]];
		d = Dimensions[mM][[1]];
		max = Table[Max[Ceiling[((1+2supp[[j]])Transpose[mM].#)[[j]] &/@{{-1/2,1/2},{1/2,1/2},{1/2,-1/2},{-1/2,-1/2}}]]+1,{j,1,d}];
		origin = max+1;
		{t2,\[CurlyPhi]MBS} = AbsoluteTiming[computeBracketSums[ck\[CurlyPhi]M,origin,mM, Compute -> "Bracket", MPAWL`Validate -> False]];
		If[StringCount[db,"Time"]>0,Print["Computing the Bracket Sums took ", t2, " seconds."]];
		saveCoefficients[Join[mM,supp],bsFile,\[CurlyPhi]MBS];
	, (* Successful: *)
		If[StringCount[db,"Time"]>0,Print["Loading Bracket Sums took ",t1," seconds."]];	
	];
	Return[{ck\[CurlyPhi]M,\[CurlyPhi]MBS}];
];


(* BS but not two files or None *)
localdlVP[g_,mM_,db_,supp_,True,orth_,s_,False] := Module[{},
	Message[delaValleePoussinMean::WrongFiles,s]; Return[$Failed]
];


(* any other case is not valid, especially when supp is not valid *)
localdlVP[g_,mM_,db_,supp_,bs_,orth_,s_,v_] := $Failed/; (
	(Length[supp]!=Dimensions[mM][[1]]) || (!VectorQ[supp,NumberQ]));


delaValleePoussinSubspaces[g_,mM_,mJ_,opts:OptionsPattern[]] :=
	localdlVPSub[g,mM,mJ,OptionValue[MPAWL`Debug], OptionValue[Orthonormalize],OptionValue[File],
	OptionValue[MPAWL`Validate]];


localdlVPSub[g_,mM_,mJ_,db_,orth_,file_,True] := 
	If[!isMatrixValid[mM] || !isMatrixValid[mJ] || (Abs[Det[mJ]]!=2) || (!isMatrixValid[Inverse[mJ].mM]),
		Return[$Failed],
		localdlVPSub[g,mM,mJ,db,orth,file,False]
	]


(* If g is a number take the pyramidFunction and the old support *)
localdlVPSub[g_,mM_,mJ_,db_,orth_,file_,False] :=
	localdlVPSub[pyramidFunction[ConstantArray[g,Dimensions[mM][[1]]],#]&,mM,mJ,db,orth,file,False] /; (NumberQ[g]);


(* Simple case: No file, no orth *)
localdlVPSub[g_,mM_,mJ_,db_,False, None,False] :=
Module[{d,dN,epsilon,mN,InvNt,adN,hM,NTg,\[Lambda]g,coeffS,coeffW,t1,BnSum,dM,InvNy,\[Epsilon],k,e},
		mN = Inverse[mJ].mM;
		d = Dimensions[mM][[1]];
		dM = patternDimension[mM, validateMatrix -> False];
		epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-dM+1;;d]]; (*Cycle lengths*)
		adN = Abs[Det[mN]];
		InvNt = Inverse[Transpose[mN]];
		hM = generatingSetBasis[Transpose[mM], Target -> "Symmetric", validateMatrix -> False];
		NTg=Transpose[mN].Complement[generatingSet[Transpose[mJ], Target -> "Symmetric", validateMatrix -> False],{{0,0}}][[1]];InvNy= Inverse[mN].((Complement[pattern[getPatternNormalform[mJ], Target -> "Symmetric"],{{0,0}}])[[1]]);
		\[Lambda]g = generatingSetBasisDecomp[NTg,Transpose[mM], Target -> "Symmetric", validateMatrix -> False];
		InvNy= Inverse[mN].((Complement[pattern[
					getPatternNormalform[mJ, validateMatrix -> False]
					,Target -> "Symmetric", validateMatrix -> False]
					,{{0,0}}])[[1]]);
		(* for these kernels, provided supp g \[SubsetEqual] [-1,1]^d this summation is enough *)
		(* BnSum[{x_,y_}] := Sum[g[{x,y}+Transpose[mJ].z]
			,{z,Flatten[Table[Table[Subscript[e,j],{j,1,d}],Evaluate[Sequence@@Table[{Subscript[e,j],-2,2},{j,1,d}]]],1]}];*)
		BnSum[x_] := Sum[g[x+Transpose[mJ].If[NumberQ[z],{z},z]]
			,{z,Flatten[Table[Table[Subscript[e,j],{j,1,d}],Evaluate[Sequence@@Table[{Subscript[e,j],-2,2},{j,1,d}]]],1]}];
		If[StringCount[db,"Text"]>0,Print["Computing the scaling function coefficients..."]];
		coeffS = ConstantArray[0,epsilon];
		t1 = AbsoluteTiming[
			Do[
				coeffS[[Sequence@@(Table[Subscript[\[Epsilon], k],{k,1,Length[epsilon]}]+1) ]] = Abs[Det[mJ]]*BnSum[InvNt.(modM[Table[Subscript[\[Epsilon], k],{k,1,Length[epsilon]}].hM,Transpose[mM],Target -> "Symmetric", validateMatrix -> False])];
				,Evaluate[Sequence@@Table[{Subscript[\[Epsilon], k],0,epsilon[[k]]-1},{k,1,Length[epsilon]}]]];
			][[1]];
		If[StringCount[db,"Time"]>0,Print["Scaling function coefficients computed in ",t1," seconds."]];
		(**)
		If[StringCount[db,"Text"]>0,Print["Computing the wavelet function coefficients..."]];
		coeffW = ConstantArray[0,epsilon];
		t1 = AbsoluteTiming[
			Do[
				coeffW[[Sequence\[NonBreakingSpace]@@\[NonBreakingSpace](Table[Subscript[\[Epsilon], k],{k,1,Length[epsilon]}]+1)]] =
					(*Abs[Det[mJ]]*)
					Exp[-2 Pi I (InvNy.(Table[Subscript[\[Epsilon], k],{k,1,Length[epsilon]}].hM))]*
					coeffS[[Sequence @@(modM[Table[Subscript[\[Epsilon], k],{k,1,Length[epsilon]}]+\[Lambda]g,DiagonalMatrix[epsilon], validateMatrix -> False]+1)]];
				,Evaluate[Sequence@@Table[{Subscript[\[Epsilon], k],0,epsilon[[k]]-1},{k,1,Length[epsilon]}]]];
			][[1]];
		If[StringCount[db,"Time"]>0,Print["Wavelet function coefficients computed in ",t1," seconds."]];
		Return[{coeffS,coeffW}];
	];


	(* Apply Orthonormalization to both results of the previous case*)
	localdlVPSub[g_,mM_,mJ_,db_,True, None,False] := 
		orthonormalizeTranslatesInSpace[#,mM,mJ,MPAWL`Debug -> db, MPAWL`Validate -> False] &/@ localdlVPSub[g,mM,mJ,db,False,None,False]


	(* Two File Strings are given, load both or compute else *)
	localdlVPSub[g_,mM_,mJ_,db_,orth_,{sFile_String,wFile_String},False] := 
	Module[{result},
		(* Try to load both *)
		If[StringCount[db,"Text"]>0, Print["Loading subspace coefficients from \[OpenCurlyDoubleQuote]",sFile,"\[CloseCurlyDoubleQuote] and \[OpenCurlyDoubleQuote]",wFile,"\[CloseCurlyDoubleQuote]..."]];
		result = {loadCoefficients[Join[mM,mJ], sFile], loadCoefficients[Join[mM,mJ], wFile]};
		If[(result[[1]]==$Failed) || (result[[2]]==$Failed),
			If[StringCount[db,"Text"]>0, Print["failed. Compute from scratch."]];
			result = localdlVPSub[g,mM,mJ,db,orth,None,False];
			(* Save *)
			saveCoefficients[Join[mM,mJ],sFile,result[[1]]];
			saveCoefficients[Join[mM,mJ],wFile,result[[2]]];
		,
			If[StringCount[db,"Text"]>0, Print["Subspace coefficients succesfully loaded."]];
		];
		Return[result];
	];


	(* All other file types *)
	localdlVPSub[g_,mM_,mJ_,db_,orth_,file_,False] := 	
	Module[{},	
		Message[delaValleePoussinSubspaces::WrongFiles, file];
		Return[$Failed];
	];


	(* All other Validate-Stuff *)
	localdlVPSub[g_,mM_,mJ_,db_,orth_,file_,v_] := $Failed;


(* ::Subsection:: *)
(*Dirichlet Kernels and subspaces*)


DirichletKernel[mM_,opts:OptionsPattern[{delaValleePoussinMean}]] :=
	delaValleePoussinMean[0,mM,opts];


DirichletKernelSubspaces[mM_,mJ_,opts:OptionsPattern[{delaValleePoussinSubspaces}]] :=
	delaValleePoussinSubspaces[0,mM,mJ,opts];


(* ::Subsection:: *)
(*The Wavelet Decomposition*)


(* g is one element (function or t), mJSet is a Vector of Matrix-Char-Vectors *)
decomposeData2D[g_, mJSetVec_, mM_, data_, opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
	decomposeData2D[ConstantArray[g, Length[mJSetVec]+1], mJSetVec, mM, data, opts] /; ((Depth[g] == 1) && VectorQ[mJSetVec, 
	 VectorQ[#, (StringQ[##] && MemberQ[localdlVPMatChars, ##]) &] &]);


(* g is a vector, mJSet is a Vector of Matrix-Chars *)
decomposeData2D[gVec_, mJSet_, mM_, data_, opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
	decomposeData2D[gVec,ConstantArray[mJSet, Length[gVec]-1], mM, data, opts] /; ((Depth[gVec] == 2) && (VectorQ[mJSet, (StringQ[#] && MemberQ[localdlVPMatChars, #]) &]));


(* level specifies the decomposition depth, g is one element and mJSet is a Vector of Matrix-Chars *)
decomposeData2D[l_,g_, mJSet_, mM_, data_, opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
	decomposeData2D[ConstantArray[g,l+1], ConstantArray[mJSet,l],mM, data, opts] /; ((Depth[g] == 1) && (VectorQ[mJSet, (StringQ[#] && MemberQ[localdlVPMatChars, #]) &]));


(* all other level specifications are not valid *)
decomposeData2D[l_,g_, mJSet_, mM_, data_, opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
	Return[$Failed] /;(IntegerQ[l] && l> 0)


(* Both are vectors of correct format *)
decomposeData2D[gVec_, mJSetVec_, mM_, data_, opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
	localDecomp[gVec,mJSetVec, mM, data, None,
		OptionValue[MPAWL`Debug], {OptionValue[DataPrefix],OptionValue[ImagePrefix], OptionValue[ImageSuffix],""},
		OptionValue[MPAWL`Validate], opts] /; ((Depth[gVec] == 2) && VectorQ[mJSetVec, 
	 VectorQ[#, (StringQ[##] && MemberQ[localdlVPMatChars, ##]) &] &] && ((Length[gVec]-1)==Length[mJSetVec]));


(* Validate *)
localDecomp[gVec_,mJSetVec_, mM_, data_, ck\[CurlyPhi]M_,db_,s_,True,
	opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
Module[{d, epsilon},
	If[!isMatrixValid[mM], Return[$Failed]];
	d = Dimensions[mM][[1]];
	If[d!= 2, Message[decomposeData2D::wrongDimensionsM, MatrixForm[mM]]; Return[$Failed]];
	epsilon = Diagonal[IntegerSmithForm[mM, ExtendedForm-> False]][[d-patternDimension[mM, validateMatrix -> False]+1;;d]];
	If[Dimensions[data]!=epsilon, Message[decomposeData2D::wrongDimensionsData, MatrixForm[epsilon], MatrixForm[Dimensions[data]]];
		Return[$Failed];
	]
	Return[localDecomp[gVec,mJSetVec,mM,data,ck\[CurlyPhi]M,db,s,False,opts]];
];


(*If we have one file suffix string, make a list of one element of that*)


localDecomp[gVec_,mJSetVec_, mM_, data_, pck\[CurlyPhi]M_,db_,{dataPre_String,imagePre_String,imageSuffix_String,path_String},False,
	opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
localDecomp[gVec,mJSetVec, mM, data, pck\[CurlyPhi]M,db,{dataPre,imagePre,{imageSuffix},path},False, opts];


localDecomp[gVec_,mJSetVec_, mM_, data_, pck\[CurlyPhi]M_,db_,{dataPre_String,imagePre_String,imageSuffixes_,path_String},False,
	opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] :=
Module[{thismJSet, mMg, mNg,mN,ck\[CurlyPhi]N,ck\[Psi]N,dataS,hatS,ScalN,dataW,hatW,wavN,ck\[CurlyPhi]M,origin,eOpts},
	eOpts = {opts}~Join~Options[decomposeData2D];
	thismJSet = First[mJSetVec]; mMg = First[gVec]; mNg = Take[gVec,{2}][[1]];
	ck\[CurlyPhi]M = pck\[CurlyPhi]M;
	If[ck\[CurlyPhi]M==None,
		ck\[CurlyPhi]M = delaValleePoussinMean[mMg,mM, MPAWL`Debug -> db, MPAWL`validateMatrix -> False,
			File -> dataPre<>path<>"ckS.dat", FilterRules[{opts}, Options[delaValleePoussinMean]]];
	];
	Do[
		mN = Inverse[dilationMatrix2D[letter]].mM;
		If[!Quiet[isMatrixValid[mN]],
			If[StringCount[db,"Text"]>0,Print["Decompose ",MatrixForm[mM]," ",If[StringLength[path]>0,"("<>path<>") ",""],"with ",letter," but ",MatrixForm[mN]," is noninteger; ending this branch."]];
		,
			If[StringCount[db,"Text"]>0,Print["Decompose ",MatrixForm[mM]," ",If[StringLength[path]>0,"("<>path<>") ",""]," with ",letter]];
			{hatS,hatW} = delaValleePoussinSubspaces[mNg,mM,dilationMatrix2D[letter],
				MPAWL`Debug -> db, MPAWL`Validate -> False, File -> {dataPre<>path<>letter<>"-S.dat",dataPre<>path<>letter<>"-W.dat"},
				FilterRules[eOpts, Options[delaValleePoussinSubspaces]]];
			origin = (Dimensions[ck\[CurlyPhi]M]+1)/2;
			{dataS,dataW} = WaveletTransformTorus[mM,dilationMatrix2D[letter],data,hatS,hatW,
					MPAWL`Debug -> db, MPAWL`Validate -> False, FilterRules[eOpts,Options[WaveletTransformTorus]]];
			origin = (Dimensions[ck\[CurlyPhi]M]+1)/2;
			ck\[Psi]N = getFourierFromSpace[hatW,ck\[CurlyPhi]M,origin,mM, MPAWL`Validate -> False];
			ck\[CurlyPhi]N = getFourierFromSpace[hatS,ck\[CurlyPhi]M,origin,mM, MPAWL`Validate -> False];
			If[OptionValue[computeWavelet],
				wavN = discretePlotFourierSeries[ConstantArray[OptionValue[PlotResolution],2],
					Sqrt[Abs[Det[mN]]]*getFourierFromSpace[dataW,ck\[Psi]N,origin,mN],origin,
					FilterRules[eOpts, Join @@ ( Options[#] & /@ {discretePlotFourierSeries,Plot,BarLegend,createBarLegend} )]];
				If[StringCount[db,"Image"]>0,Print[wavN]];
				Do[
				If[StringCount[db,"Text"]>0,Print["Exporting to File \[OpenCurlyDoubleQuote]",imagePre<>"Wavelet-"<>path<>letter<>imageSuf,"\[CloseCurlyDoubleQuote]"]];			
				Export[imagePre<>"Wavelet-"<>path<>letter<>imageSuf,
					 Show[wavN,Sequence@@FilterRules[{opts}, Options[Show]~Join~Options[Graphics]]],
						Sequence@@FilterRules[eOpts,FilterRules[Options[Export]~Join~Options[Rasterize], Except[ImageSize]]]];
				,{imageSuf,imageSuffixes}];
			];
			If[OptionValue[computeScale],
				ScalN = discretePlotFourierSeries[ConstantArray[OptionValue[PlotResolution],2],
						Sqrt[Abs[Det[mN]]]*getFourierFromSpace[dataS,ck\[CurlyPhi]N,origin,mN],origin,
						FilterRules[eOpts, Join @@ ( Options[#] & /@ {discretePlotFourierSeries,Plot,BarLegend,createBarLegend} )]];
				If[StringCount[db,"Image"]>0,Print[ScalN]];
				Do[
				If[StringCount[db,"Text"]>0,Print["Exporting to File \[OpenCurlyDoubleQuote]",imagePre<>"Wavelet-"<>path<>letter<>imageSuf,"\[CloseCurlyDoubleQuote]"]];			
				Export[imagePre<>"Scale-"<>path<>letter<>imageSuf,
					Show[ScalN,Sequence@@FilterRules[{opts}, Options[Show]~Join~Options[Graphics]]],
						Sequence@@FilterRules[eOpts,FilterRules[Options[Export]~Join~Options[Rasterize], Except[ImageSize]]]];
				,{imageSuf,imageSuffixes}];
			];
			If[Length[Rest[mJSetVec]]>0, (* Still at least one level to go *)
				localDecomp[Rest[gVec],Rest[mJSetVec],mN,dataS,ck\[CurlyPhi]N,db,{dataPre,imagePre,imageSuffixes,path<>letter},False,opts]
			];
		];
		,{letter, thismJSet}
	]
]/; (VectorQ[imageSuffixes,StringQ]);


(* all other cases, i.e. path not a String or such things*)
localDecomp[gVec_,mJSetVec_, mM_, data_, ck\[CurlyPhi]M_,db_,s_,v_,
	opts:OptionsPattern[{decomposeData2D,Export,Show,Plot,BarLegend,createBarLegend,discretePlotFourierSeries}]] := $Failed;


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
