(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Visualization*)


(* ::Subsubtitle:: *)
(*Functions that provide several visualization possibibilies for Fourier Series and pattern / generating Set points*)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		25.03.2011*)
(*Last Changed: 	03.03.2013*)


BeginPackage["MPAWL`Visualization`",
{
"MPAWL`Basics`",
"MPAWL`Pattern`",
"MPAWL`Transforms`"
}];


(* ::Section:: *)
(*Global Function Declaration*)


(* ::Text:: *)
(*Declaration of all functions that this package provides for global usage*)


(* ::Subsection:: *)
(*Functions for Visualization of Fourier Sums.*)


discretePlotFourierSeries::usage="discretePlotFourierSeries[resolution,coefficients,origin]

For an array of Fourier coefficients, where origin represents the index of the
origin, an Image of size resolution is produced, showing the function either
in color (where red indicated positive and blue negative values), grayscale or
absolute value.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

ReturnVal \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]AbsoluteImage\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]ColorImage\[CloseCurlyDoubleQuote]
	Specify, which kind of Image should be produced
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.
ColorLegend \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | True
	display a range of values and their corresponding colors to the right of
	the plot.

The Image is included in a plot, and hence accepts all its options,
providing some non-standard Options as new standards. The same holds for the
createBarLegend and hence BarLegend the ColorLegend option consists of.";


discretePlotFourierSeriesDiff::usage="discretePlotFourierSeriesDiff(resolution,coefficients,origin,f,(options)]

Perform the same production of an image as in discretePlotFourierSeries[], but
additionaly subtract the values of the function f at each pixel obtained.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)
ReturnVal \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]AbsoluteImage\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]ColorImage\[CloseCurlyDoubleQuote]
	Specify, which kind of Image should be produced
Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.
ColorLegend \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | True
	display a range of values and their corresponding colors to the right of
	the plot.

The Image is included in a plot, and hence accepts all its options,
providing some non-standard Options as new standards. The same holds for the
BarLegend the ColorLegend consists of.";


discretePlotFourierSeries::wrongResolution = "The resolution `1` does not match the Dimensions of the coefficients `2`";
discretePlotFourierSeries::wrongOrigin = "The origin `1` does not match the Dimensions of the coefficients `2`";


Options[discretePlotFourierSeries] = {
MPAWL`Debug -> "None", ReturnVal -> "AbsoluteImage",ColorLegend -> False,
(* to override Plot *)
Axes -> False, Frame -> False, AxesOrigin -> {-\[Pi],-\[Pi]},AxesLabel -> {"x","y"},
TicksStyle -> Directive[8,Plain],
FrameTicks -> {{Range[-\[Pi],\[Pi],\[Pi]/2],None},{Range[-\[Pi],\[Pi],\[Pi]/2],None}},
PlotRangePadding -> 0,AspectRatio -> 1};


(* ::Subsection:: *)
(*Plot on Torus*)


plotOnTorus::usage = "PlotPatOnTorus[in]

Plots the input on a torus, which can be used for various input types:
	1) a Matrix, that is integral and has full rank, the its pattern is plotted
	2) an arbitrary set of points on the 2\[Pi]-periodic torus
	3) a function, that is assumed to be defined at least on [0,2\[Pi])^2

In order to change the color of the surface, change ColorFunction.
In order to change the color of the plottet points (cases 1 and 2),
change the PLotStyle

Further PointSize can be specified to raise the points above the torus.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote])
	to produce intermediate results, indicate progress and display computation
	times.";


plotOnTorus::noSetofPoints = "The input seems to be neither a set of points nor a set of sets of points";


Options[plotOnTorus] = {Boxed -> False,
Axes -> False,
PlotRange -> All,
Background -> None,
ColorFunction -> Function[{x,y,z},RGBColor[1,1,1,.83]],
MaxRecursion -> 5,
Mesh -> {3,7},
MeshStyle -> RGBColor[.5,.5,.5,.75],
BoundaryStyle -> RGBColor[0,75/255,90/255],
ViewPoint -> {2,0,2.5}, 
ViewVertical -> {0,0,1},
PointSize-> 0.02,
PlotStyle ->  Directive[PointSize[0.02],Opacity[.75],RGBColor[0,75/255,90/255]],
MPAWL`Debug -> "None"};


createBarLegend::usage="createBarLegend[min,max]

Create a BarLegend ranging from min to max for all IMageTypes of discretePlotFourierSeries.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

ScientificNumners \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | True
	the Legend displays it's values in Scientific Form, if set to true.
AddToLegendExponent \[Rule] \!\(\*StyleBox[\"False\",\nFontSlant\[Rule]\"Italic\"]\) | natural Number
	for some cases, the scientific numbers fail, hence the data has to be rescaled.
	This option fixes the display.";


Options[createBarLegend] = 
{(*LegendMargins -> {{2,2},{-44,2}},*)
LabelStyle -> Directive[8,Plain,Black],
ScientificNumbers -> True,
AddToLegendExponent -> 0,
MPAWL`Debug -> "None",
ReturnVal -> "AbsoluteImage"};


(* ::Section:: *)
(*Function Definition*)


(* ::Text:: *)
(*All function definitions are inside the private part, so only those functions are globally visible and thus usable, that are declared above via ::usage*)
(*And even their Error Messages are not visible to public*)


Begin["`Private`"];


localPaddAndFourierTransform::usage = "localPaddAndFourierTransform[resolution, coefficients, origin]
Expand the array of coefficients to coincide with resolution and perform a
Fourier transform on that array.

\!\(\*StyleBox[\"Options\",FontWeight\[Rule]\"Bold\"]\)

Debug \[Rule] \!\(\*StyleBox[\"\[OpenCurlyDoubleQuote]None\[CloseCurlyDoubleQuote]\",\nFontSlant\[Rule]\"Italic\"]\) | \[OpenCurlyDoubleQuote]Text\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Time\[CloseCurlyDoubleQuote] | \[OpenCurlyDoubleQuote]Image\[CloseCurlyDoubleQuote]
	or any combination of these Words in one String (i.e. concatenated via \[OpenCurlyDoubleQuote]&\[CloseCurlyDoubleQuote]
	to produce intermediate results, indicate progress and display computation
	times.

Validate \[Rule] \!\(\*StyleBox[\"True\",\nFontSlant\[Rule]\"Italic\"]\) | False
		whether to perform a check (via isMatrixValid[mM]) on the matrix mM
		and the check, whether the Origin is in Range.";


Options[localPaddAndFourierTransform] = {MPAWL`Debug -> "None", MPAWL`Validate -> True}; 


Options[createImageFromArray] = {ReturnVal -> "AbsoluteImage", MPAWL`Debug -> "None"};


(* ::Subsection:: *)
(*Helpers*)


localPaddAndFourierTransform[resolution_, coefficients_,origin_,opts:OptionsPattern[]] :=
	localPAndFT[resolution, coefficients, origin, OptionValue[MPAWL`Debug], OptionValue[MPAWL`Validate]];

localPAndFT[resolution_, coefficients_, origin_, db_, True] := Module[{},
	If[(Length[resolution]!=2) || (Length[origin]!=2) 
		|| (Length[Dimensions[coefficients]] != 2),
	Message[discretePlotFourierSeries::wrongResolution,
		MatrixForm[resolution],MatrixForm[Dimensions[Dimensions[coefficients]][[1]]]];
	Return[$Failed];
	];
	Return[localPAndFT[resolution, coefficients, origin, db, False]];
];


localPAndFT[resolution_, coefficients_, origin_, db_, False] :=
Module[{k,image,i,fouriert,t},
	If[StringCount[db,"Text"]>0,Print["Padding Data..."]];
	image = ConstantArray[0,resolution];
	t = {0,0};
	Do[
		t[[1]] = k1+origin[[1]]-1;
		If[k1 > Dimensions[coefficients][[1]]- origin[[1]]+1, t[[1]] -= resolution[[1]];];
		t[[2]] = k2+origin[[2]]-1;
		If[k2 > Dimensions[coefficients][[2]]- origin[[2]]+1, t[[2]] -= resolution[[2]];];
		If[isIndexInRange[coefficients, t], image[[k1,k2]] = coefficients[[Sequence@@t]]];
	,{k1,1,resolution[[1]]},{k2,1,resolution[[2]]}];
	If[StringCount[db,"Text"]>0,Print["Fourier transforming..."]];
	fouriert = FourierTransformTorus[DiagonalMatrix[resolution],
		Sqrt[Times[Sequence@@resolution]]*image,
		MPAWL`Validate -> False];
	If[StringCount[db,"Text"]>0,Print["Shifting back..."]];
	image = ConstantArray[0,resolution];
	Do[
	image[[k1,k2]] = fouriert[[Sequence@@(Mod[{k1,k2}+Ceiling[resolution/2],resolution]+1)]];
	,{k1,1,resolution[[1]]},{k2,1,resolution[[2]]}];
	Return[image];
];


localPAndFT[resolution_, coefficients_, origin_, db_, v_] := $Failed;


createImageFromArray[image_,opts:OptionsPattern[]] := 
	localCreateImage[image,OptionValue[MPAWL`Debug], OptionValue[ReturnVal]];


localCreateImage[image_,db_,"AbsoluteImage"] := Module[{lImage},
	If[StringCount[db,"Text"]>0,Print["Creating the absolute image..."]];
	lImage = Abs[image];
	If[StringCount[db,"Text"]>0,Print["The range of values is ",ScientificForm[Min[Chop[N[Re[lImage]]]]]," to ",ScientificForm[Max[Chop[N[Re[lImage]]]]],"."];];
	(* Normalize *)
	lImage = lImage - Min[lImage];
	If[Max[lImage]!= 0,lImage = lImage / Max[lImage];];
	Return[Image[Chop[lImage]]]
]


localCreateImage[image_,db_,"Image"] := Module[{lImage},
	If[StringCount[db,"Text"]>0,Print["Creating the image..."]];
	If[StringCount[db,"Text"]>0,Print["The range of values is ",ScientificForm[Min[Chop[N[Re[image]]]]]," to ",ScientificForm[Max[Chop[N[Re[image]]]]],"."];];
	(* Normalize *)
	lImage = image - Min[image];
	If[Max[lImage]!= 0,lImage = lImage / Max[lImage];];
	Return[Image[Chop[lImage]]]
]


localCreateImage[image_,db_,"ColorImage"] :=
Module[{posImage,negImage,resolution,norm,lImage},
	If[StringCount[db,"Text"]>0,Print["Creating Color Image..."]];
	resolution = Dimensions[image];
	If[StringCount[db,"Text"]>0,Print["The range of Values is from ",ScientificForm[Min[Chop[N[Re[image]]]]]," to ",ScientificForm[Max[Chop[N[Re[image]]]]],"."];];
	posImage = Table[Max[image[[x,y]],0],{x,1,resolution[[1]]},{y,1,resolution[[2]]}];
	negImage = Table[Min[0,image[[x,y]]],{x,1,resolution[[1]]},{y,1,resolution[[2]]}];
	norm = Max[Max[posImage], Max[negImage]];
	posImage = posImage / If[norm==0,1,norm];
	negImage = negImage / If[norm==0,1,norm];(* Shift both *)
	posImage = posImage + (1-Max[posImage]);
	negImage = negImage + (1-Max[negImage]);
	(*Hard Coded Colors for Pos and Neg: Red & Blue*)
	Return[
		Image[  Table[1 - (1-negImage[[x,y]])*{1,1,0} - posImage[[x,y]]*{0,1,1},{x,1,resolution[[1]]},{y,1,resolution[[2]]}]  ]
		];
];


localCreateImage[image_,db_,img_] := $Failed;


createBarLegend[min_,max_,opts:OptionsPattern[{createBarLegend, BarLegend}]] :=
Module[{exp,formatLabel, barlegendOpts,x,y},
exp = Floor[Log[10,Max[Abs[min],Abs[max]]]];
formatLabel[x_] := x/.{NumberForm[y_,{w_,z_}] :> NumberForm[y/(10^(exp)),{1,2}]};
barlegendOpts = FilterRules[{opts}, Options[BarLegend]];
If [OptionValue[ScientificNumbers],
	barlegendOpts = FilterRules[{barlegendOpts}, Except[{LegendFunction,LegendLabel}]];
	barlegendOpts = Append[barlegendOpts,
		{LegendFunction->formatLabel,
		LegendLabel->
		If[exp+OptionValue[AddToLegendExponent]!=0, Placed[DisplayForm[SuperscriptBox[ToString["  \[Times] 10"],(exp+OptionValue[AddToLegendExponent])]],Bottom]
		, None]
			}]
];
Return[
If[(OptionValue[ReturnVal] == "ColorImage"),
	BarLegend[
		{RGBColor[{1-Max[1-2(#-min)/(max-min),0],
		1-Max[1-2(#-min)/(max-min),0]-Max[-1+2(#-min)/(max-min),0],
		1-Max[-1+2(#-min)/(max-min),0]}] &,{min,max}},
		Sequence@@barlegendOpts
			]
	,
	(* Graysclae *)
	BarLegend[{RGBColor[(1-(#-min)/(max-min)){1,1,1}]&,{min,max}},Sequence@@barlegendOpts]
]];
];


(* ::Subsection:: *)
(*Visualization of Fourier Sums*)


discretePlotFourierSeriesDiff[resolution_, coefficients_,origin_,f_,opts:OptionsPattern[{discretePlotFourierSeries,createBarLegend, BarLegend}]]:= 
Module[{image,dim,i,min,max,minmax,j,t1,dataImg,remainingRulles,k,j},
dim = Length[resolution];
t1 = AbsoluteTiming[
	image = localPaddAndFourierTransform[resolution,coefficients,origin,MPAWL`Debug-> OptionValue[MPAWL`Debug], MPAWL`Validate-> False];
	Do[
		image[[Sequence@@Table[k[j],{j,1,dim}]\[NonBreakingSpace]]] = 
		f[2\[Pi] (Table[(k[j] - resolution[[j]]/2)/resolution[[j]],{j,1,dim}]) ] - 
		image[[Sequence@@Table[k[j],{j,1,dim}]\[NonBreakingSpace]]];
		,Evaluate[Sequence@@ Table[{k[l],1,resolution[[l]]},{l,1,dim}]]
	];
	image = Chop[N[Re[image]]];
	min = Min[image];max=Max[image];
	If[OptionValue[ReturnVal] == "AbsoluteImage",min = Min[Abs[image]]; max = Max[Abs[image]]];
	minmax = Max[Abs[min],Abs[max]];
	If[OptionValue[ReturnVal] == "ColorImage",min = -minmax; max = minmax];
	image = createImageFromArray[image,FilterRules[{opts},Options[createImageFromArray]]];
][[1]];
If[StringCount[OptionValue[MPAWL`Debug],"Time"]>0,Print["Creating the image took ",t1," seconds"]];
If[OptionValue[ReturnVal] == "Matrix",Return[image]];
dataImg = ArrayPlot[
If[OptionValue[ReturnVal]=="ColorImage", Map[RGBColor,ImageData[ImageRotate[image]],{2}],ImageData[ImageRotate[image]]],
DataRange-> {{-\[Pi],\[Pi]},{-\[Pi],\[Pi]}},PlotLegends-> If[OptionValue[ColorLegend],
	Placed[createBarLegend[min,max,	FilterRules[{opts},Options[createBarLegend]~Join~Options[BarLegend]]],{1.01,.48}],None]];
	remainingRules = FilterRules[Options[discretePlotFourierSeries],Except[{opts}]];
	Return[Show[Plot[0,{x,-\[Pi],\[Pi]},PlotRange->{-\[Pi],\[Pi]},PlotStyle-> Directive[Opacity[0]],
					Evaluate[Sequence@@FilterRules[remainingRules~Join~{opts},Options[Plot]]]
					],dataImg]];
];


discretePlotFourierSeries[resolution_, coefficients_,origin_,opts:OptionsPattern[{discretePlotFourierSeries,createBarLegend, BarLegend,Plot}]]:= 
Module[{image,dim,i,min,max,minmax,j,t1,dataImg, remainingRules},
	dim = Length[resolution];
	t1 = Timing[
		image = localPaddAndFourierTransform[resolution,coefficients,origin,FilterRules[{opts},Options[localPaddAndFourierTransform]]];
		image = Chop[N[Re[image]]];
		min = Min[image];max=Max[image];
		If[OptionValue[ReturnVal] == "AbsoluteImage",min = Min[Abs[image]];
		max = Max[Abs[image]]];
		minmax = Max[Abs[min],Abs[max]];
		If[OptionValue[ReturnVal] == "ColorImage",min = -minmax; max = minmax];
		image = createImageFromArray[image,FilterRules[{opts},Options[createImageFromArray]]];
	][[1]];
If[OptionValue[ReturnVal] == "Matrix",Return[image]];
dataImg = ArrayPlot[
If[OptionValue[ReturnVal]=="ColorImage", Map[RGBColor,ImageData[ImageRotate[image]],{2}],ImageData[ImageRotate[image]]],
DataRange-> {{-\[Pi],\[Pi]},{-\[Pi],\[Pi]}},PlotLegends-> If[OptionValue[ColorLegend],
	Placed[createBarLegend[min,max,	FilterRules[{opts},Options[createBarLegend]~Join~Options[BarLegend]]],{1.03,.5}],None]];
	remainingRules = FilterRules[Options[discretePlotFourierSeries],Except[{opts}]];
	Return[Show[Plot[0,{x,-\[Pi],\[Pi]},PlotRange->{-\[Pi],\[Pi]},PlotStyle-> Directive[Opacity[0]],
					Evaluate[Sequence@@FilterRules[remainingRules~Join~{opts},Options[Plot]]]
					],dataImg]];
];


(* ::Subsection:: *)
(*Plot on Torus*)


(* Private taken from
http://stackoverflow.com/questions/3736942/test-if-an-expression-is-a-function
checks whether the input is function *)
FunctionQ[_Function | _InterpolatingFunction | _CompiledFunction] = True;
FunctionQ[f_Symbol] := Or[
  DownValues[f] =!= {}, 
  MemberQ[ Attributes[f], NumericFunction ]]
FunctionQ[_] = False;


plotOnTorus[in_,opts:OptionsPattern[{plotOnTorus,ParametricPlot3D,ListPointPlot3D}]] :=
Module[{bg,valplot,pts,torusf,x,y},
	(*Set Standardvalues for ListPointPlot3D *)
	SetOptions[ListPointPlot3D,Evaluate[Sequence@@FilterRules[Options[plotOnTorus],Complement[Options[ListPointPlot3D],Options[ListPointPlot3D,ColorFunction]]]]];
	(*Set Standardvalues for ParametricPlot3D despite Colorfunction*)
	SetOptions[ParametricPlot3D,Evaluate[Sequence@@FilterRules[Options[plotOnTorus],Options[ParametricPlot3D]]]];
	If[Quiet[isMatrixValid[in]],
		pts = pattern[getPatternNormalform[in], Target -> "Symmetric"];
		If[StringCount[OptionValue[MPAWL`Debug],"Text"]>0,Print["Interpreting the input as a Matrix"]];
		valplot = ListPointPlot3D[
		Table[{Cos[k[[2]]2Pi] (2+(1+OptionValue[PointSize])Cos[k[[1]]2Pi]),Sin[k[[2]]2Pi] (2+(1+OptionValue[PointSize])Cos[k[[1]]2Pi]),(1+OptionValue[PointSize])Sin[k[[1]]2Pi] },{k,pts}]
		,	Evaluate[Sequence@@Complement[(*Exclude Colorfunction *)FilterRules[{opts},Options[ListPointPlot3D]],FilterRules[{opts},{ColorFunction}]]]
			];
		bg = ParametricPlot3D[{Cos[t] (2+Cos[u]),Sin[t] (2+Cos[u]),Sin[u]}, {t,-\[Pi],\[Pi]},{u,-\[Pi],\[Pi]},
			Evaluate[Sequence@@FilterRules[{opts},Options[ParametricPlot3D]]]];
	Return[Show[bg,valplot,PlotRange->Automatic]];
		,
		If[FunctionQ[in],
			If[StringCount[OptionValue[MPAWL`Debug],"Text"]>0,Print["Interpreting the input as a function."]];
			torusf[k_] :=
				{
					Cos[k[[2]]] (2+(1+in[k])Cos[k[[1]]]),
					Sin[k[[2]]] (2+(1+in[k])Cos[k[[1]]]),
					(1+in[k])*Sin[k[[1]]]
				};
			Return[ParametricPlot3D[torusf[{x,y}],{x,-\[Pi],\[Pi]},{y,-\[Pi],\[Pi]},
Evaluate[Sequence@@FilterRules[{opts},Options[ParametricPlot3D]]]]];
			,
			If[Length[Dimensions[in]]==2,
				If[StringCount[OptionValue[MPAWL`Debug],"Text"]>0,Print["Interpreting the input as a set of points"]];
				pts = in;
		valplot = ListPointPlot3D[
		Table[{Cos[k[[2]]] (2+(1+OptionValue[PointSize])Cos[k[[1]]]),Sin[k[[2]]] (2+(1+OptionValue[PointSize]/2)Cos[k[[1]]]),(1+OptionValue[PointSize])Sin[k[[1]]] },{k,pts}],
			Evaluate[Sequence@@Complement[(*Exclude Colorfunction *)FilterRules[{opts},Options[ListPointPlot3D]],FilterRules[{opts},{ColorFunction}]]]
			];
			, (* not 2 *)
			If[Length[Dimensions[in]]==3,
				If[StringCount[OptionValue[MPAWL`Debug],"Text"]>0,Print["Interpreting the input as a set of sets of points"]];
				pts = in;
			valplot = ListPointPlot3D[
			Table[
				Table[{Cos[k[[2]]] (2+(1+OptionValue[PointSize])Cos[k[[1]]]),Sin[k[[2]]] (2+(1+OptionValue[PointSize]/2)Cos[k[[1]]]),(1+OptionValue[PointSize])Sin[k[[1]]] },{k,onepts}]
			,{onepts,pts}],
			Evaluate[Sequence@@Complement[(*Exclude Colorfunction *)FilterRules[{opts},Options[ListPointPlot3D]],FilterRules[{opts},{ColorFunction}]]]
			];
			,
			Message[plotOnTorus::noSetofPoints];Return[$Failed];
			]];
		bg = ParametricPlot3D[{Cos[t] (2+Cos[u]),Sin[t] (2+Cos[u]),Sin[u]}, {t,-\[Pi],\[Pi]},{u,-\[Pi],\[Pi]},
		Evaluate[Sequence@@FilterRules[{opts},Options[ParametricPlot3D]]]];
		Return[Show[bg,valplot,PlotRange->Automatic]];
		];
	];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ]


EndPackage[ ]
