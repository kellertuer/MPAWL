(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Box Splines*)


(* ::Subsubtitle:: *)
(*Adaption of an algorithm proposed by L. Kobbelt (1997) from MatLab to Mathematica*)
(*  *)


(* ::Text:: *)
(*This part of the Library should not be included by itself. Instead the whole Library should be loaded  by*)
(*using Needs["MPAWL`"].*)


(* ::Text:: *)
(*This file contains an implementation of the algorithm proposed by L. Kobbelt for a stable evaluation of Box Splines,*)
(*   though still working with the normal [0, 1)^d characteristic function instead of the formalism used*)
(*to propose properties of the wavelets.*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		13.11.2012*)
(*Last Changed: 	02.03.2013*)


BeginPackage["MPAWL`BoxSplines`"];


(* ::Section:: *)
(*Global Function Declaration*)


evaluateBoxSpline::usage="evaluateBoxSpline[\[CapitalXi],x]
evaluateBoxSpline[\[CapitalXi],\[Nu],x]
evaluateBoxSpline[x]

Evaluate s Box-Spline at the point or set of points x, each d-dimensional.

The first \[CapitalXi] defines the Box-Spline using a multiset of vectors.
The second by set of vectors \[CapitalXi] and its multiplicities \[Nu].
If there was already a call of this function, the last one uses the last
initialized Box-Spline for evaluation.";


evaluateBoxSpline::PointUnsuitable = "The Point `1` you specified lies not inside the definition Range of the Box-Spline (which is the \!\(\*SuperscriptBox[\(\[DoubleStruckCapitalR]\), \(`2`\)]\)).";
evaluateBoxSpline::NotInitialized = "Only specifying a Point is not enough if there is no initialized Box Spline.";


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


(* ::Text:: *)
(*These internal values get initialized when the package is read and store the last evaluated BoxSpline in intermediate Values*)


k=0;  (* Number of rows in X *)
U={}; (* Hashtable N in the code of Kobbelt *)
u={};
d=0;  (* s in MATLAB *)
X={}; (* distinct Vecotrs in \[CapitalXi] *)
nu={}; (*multiplicities of vectors in \[CapitalXi] that are uniquely stored in X *)
initialized = False;
Y={};
p={};


localInitBoxSpline::usage = "If the set os vectors is given in the form of \[CapitalXi] without multiplicities,
this method performs the initialization and returns the success of initialization";
localInitBoxSpline::wrongDirections = "The BoxSpline Vectors `1`are not in correct shape";
localInitBoxSpline[\[CapitalXi]_] := Module[{lX,lnu},
lX = Union[\[CapitalXi]];
lnu = ConstantArray[0,Dimensions[lX][[1]]];
Do[
	lnu[[i]] = Count[\[CapitalXi],lX[[i]]];
	,{i,1,Dimensions[lX][[1]]}
];
localInitBoxSpline[lX,lnu];
]
localInitBoxSpline[\[CapitalXi]_,\[Nu]_] := Module[{},
X = \[CapitalXi]; nu = \[Nu];
If [(Dimensions[Dimensions[X]] != 2),
	Message[localInitBoxSpline::wrongDirections,X];
	initialized=False;
	Return[False];
];
{k,d} = Dimensions[X];
u = Table[2^j,{j,0,k-1}];
Y = Extract[X,Position[nu,Except[0,_Integer],1]]; (* Extract all at least once occuring Elements from X *) 
(* BackslashOperator is nicer in Matlab, because it automatically solves the LGS for each column of Y*)
Y = Table[LinearSolve[Transpose[Y].Y,Y[[j]]],{j,1,Dimensions[Y][[1]]}]; 
U = localBoxSplineNorm[d-1,k,ConstantArray[0,k]];
initialized=True;
Return[True];
];


evaluateBoxSpline[\[CapitalXi]_,x_] := Module[{},
If [!localInitBoxSpline[\[CapitalXi]],
	Return[$Failed];
];
If [d!= Dimensions[x][[1]],
	Message[evaluateBoxSpline::PointUnsuitable,x,d];
	Return[$Failed];
];
p = x;
Return[BoxSplineRecursion[nu, ConstantArray[0,k],Y,Y.x]]
]/; (Length[Dimensions[x]]==1)


evaluateBoxSpline[\[CapitalXi]_,p\[Nu]_,x_] := Module[{},
If [!localInitBoxSpline[\[CapitalXi],p\[Nu]],
	Return[$Failed];
];
p = x;
Return[BoxSplineRecursion[nu, ConstantArray[0,k],Y,Y.x]]
]/;(Length[Dimensions[x]]==1)


evaluateBoxSpline[x_] := Module[{},
	If[!initialized,
		Message[evaluateBoxSpline::NotInitialized];
		Return[$Failed]
	];
	p = x;
	(*Else all values nu, X, Y and so on exist, start recursive evaluation*)
	Return[BoxSplineRecursion[nu, ConstantArray[0,k],Y,Y.x]]
]/;(Length[Dimensions[x]]==1)


evaluateBoxSpline[\[CapitalXi]_,x_] := Module[{},
If [!localInitBoxSpline[\[CapitalXi]],
	Return[$Failed];
];
If [d!= Dimensions[x][[2]],
	Message[evaluateBoxSpline::PointUnsuitable,x,d];
	Return[$Failed];
];
p = x;
Return[BoxSplineRecursionMult[nu, ConstantArray[0,k],Y,Y.#&/@x]]
]/;(Length[Dimensions[x]]==2) (*set of points*)


evaluateBoxSpline[\[CapitalXi]_,p\[Nu]_,x_] := Module[{},
If [!localInitBoxSpline[\[CapitalXi],p\[Nu]],
	Return[$Failed];
];
p = x;
Return[BoxSplineRecursionMult[nu, ConstantArray[0,k],Y,Y.#&/@x]]
]/;(Length[Dimensions[x]]==2) (*set of points*)


evaluateBoxSpline[x_] := Module[{},
	If[!initialized,
		Message[evaluateBoxSpline::NotInitialized];
		Return[$Failed]
	];
p = x;
Return[BoxSplineRecursionMult[nu, ConstantArray[0,k],Y,Y.#&/@x]]
]/;(Length[Dimensions[x]]==2) (*set of points*)


localBoxSplineNorm::usage = "localBoxSplineNorm[t,k,M]

Precomputation of Normalvectors in the Package-Global N, based on L. Knobbelt,
where 
	t ist the number og rows selected before
	k is the next row of X to be considered
	M is the Bitvector indicating the selected rows of X
";


localBoxSplineNorm[t_,pk_,M_] := Module[{nspace},
	U = ConstantArray[0,{d,2^pk}];
If [pk >= t,
	If [t > 0,
		U = Join[localBoxSplineNorm[t,pk-1,M],localBoxSplineNorm[t-1,pk-1,M+UnitVector[k,pk]],2];
	,
	(*Normal vector is orthogonal to all selected rows ... and unique*)
			nspace = NullSpace[Extract[X,Position[M,Except[0,_Integer],1]]][[1]];
			U[[All,1]] = nspace/Norm[nspace];
		];
	];
Return[U];
];


BoxSplineRecursion::usage = "BoxSplineRecursion[n,m,Y,t]

Recursive computation based on L. Knobbelt (1997), where

	n holds the number of multiplicities of the vectors in X,
	m is the current position in recursion tree
	Y is the Matrix to compute least norm representation of p
	t least norm representation of p";


BoxSplineRecursion[n_,m_,Y_,t_] := Module[{b,j,i,sumn,nn,mm,Z,v,z,lp,lq,NN},
(* recursion case more than d vectors left *)
	sumn = Plus@@n;
	If[ (sumn) > d, (*recursdion case *)
	b=0; j=1;
	(*Sum over remaining Elements of X *)
	For [i = 1, i <= k, i++,
		nn = n-UnitVector[k,i];
		mm = m+UnitVector[k,i];
		(* Recursive calls *)
		If [n[[i]] > 1,
			b = b + t[[j]]*BoxSplineRecursion[nn,m,Y,t];
			b = b + (n[[i]]-t[[j]])*BoxSplineRecursion[nn,mm,Y,
							t-(X[[i]].Transpose[Y])];
			j=j+1;
			, (*else ( < 1*)
			If [n[[i]] > 0,
				Z = Extract[X,Position[nn,Except[0,_Integer],1]]; (* Extract all at least once occuring Elements from X *) 
				If [MatrixRank[Z] == d,
					Z = Table[LinearSolve[Transpose[Z].Z,Z[[j]]],{j,1,Dimensions[Z][[1]]}];
					b = b + t[[j]]*BoxSplineRecursion[nn,m,Z,
									(p-(m.X)).Transpose[Z]];
					b = b + (n[[i]]-t[[j]])*BoxSplineRecursion[nn,mm,Z,
									(p-(mm.X)).Transpose[Z]];
				];
				j=j+1;
			];
		];
	];
	b = b/(sumn-d);
	Return[b];
, (*else end of recursion *)
	b=1;
	v = Position[n,Except[0,_Integer],1]; (* Extract all at least once occuring Elements from X *) 
	z = p-(m.X); (* delayed translation *)
	(* Check against all Hyperplanes*)
	For[i=1,i <= d, i++,
		(* normal vector to ith hyperplane *)
		NN = U[[All,1+u.(n-UnitVector[k,v[[i,1]]])]];
		(*lp: relevant sides of the hypercube *)
		lp = Extract[X,v[[i]]].NN; (* X[[v[[i]]]].NN;*)
		(* projection of translated point onto hyperplane*)
		lq = z.NN;
		b = Min[b,If[(((lp>0)&&(lq<0))||((lp<0)&&(lq >= 0))),0,1]];
		lq = (p-((m+UnitVector[k,v[[i,1]]]).X)).NN; (*TODO: Fix for more than on p*)
		b= Min[b,If[(((lp>0)&&(lq >= 0))||((lp<0)&&(lq<0))),0,1]];
	]; 
	(* Normalization *)
	b = b/(Abs[Det[Extract[X,v[[1;;d]]]]]);
	Return[b];
];
];


BoxSplineRecursionMult[n_,m_,Y_,t_] := Module[{b,j,i,sumn,nn,mm,Z,v,z,lp,lq,lq2,NN},
(* recursion case more than d vectors left *)
	sumn = Plus@@n;
	If[ (sumn) > d, (*recursdion case *)
	b=ConstantArray[0,Length[t]]; j=1;
	(*Sum over remaining Elements of X *)
	For [i = 1, i <= k, i++,
		nn = n-UnitVector[k,i];
		mm = m+UnitVector[k,i];
		(* Recursive calls *)
		If [n[[i]] > 1,
			b = b + t[[All,j]]*BoxSplineRecursion[nn,m,Y,t];
			b = b + (n[[i]]-t[[All,j]])*BoxSplineRecursionMult[nn,mm,Y,
							(#-(X[[i]].Transpose[Y]) & /@t)];
			j=j+1;
			, (*else ( < 1*)
			If [n[[i]] > 0,
				Z = Extract[X,Position[nn,Except[0,_Integer],1]]; (* Extract all at least once occuring Elements from X *) 
				If [MatrixRank[Z] == d,
					Z = Table[LinearSolve[Transpose[Z].Z,Z[[j]]],{j,1,Dimensions[Z][[1]]}];
					b = b + t[[All,j]]*BoxSplineRecursionMult[nn,m,Z,
									(#-(m.X)).Transpose[Z] & /@p];
					b = b + (n[[i]]-t[[All,j]])*BoxSplineRecursionMult[nn,mm,Z,
									(#-(mm.X)).Transpose[Z] & /@p];
				];
				j=j+1;
			];
		];
	];
	b = b/(sumn-d);
	Return[b];
, (*else end of recursion *)
	b=ConstantArray[1,Length[t]];
	v = Position[n,Except[0,_Integer],1]; (* Extract all at least once occuring Elements from X *) 
	z = #-(m.X) & /@p; (* delayed translation of all points*)
	(* Check against all Hyperplanes*)
	For[i=1,i <= d, i++,
		(* normal vector to ith hyperplane *)
		NN = U[[All,1+u.(n-UnitVector[k,v[[i,1]]])]];
		(*lp: relevant sides of the hypercube *)
		lp = Extract[X,v[[i]]].NN; (* X[[v[[i]]]].NN;*)
		(* projection of translated point onto hyperplane*)
		lq = #.NN & /@z;
		lq2 = (#-((m+UnitVector[k,v[[i,1]]]).X)).NN & /@p;
		b = Table[
			Min[
				Min[b[[k]],If[(((lp>0)&&(lq[[k]]<0))||((lp<0)&&(lq[[k]]>= 0))),0,1]],
			If[(((lp>0)&&(lq2[[k]] >= 0))||((lp<0)&&(lq2[[k]]<0))),0,1]
			],
			{k,1,Length[p]}
			];
	];
	(* Normalization *)
	b = b/(Abs[Det[Extract[X,v[[1;;d]]]]]);
	Return[b];
];
];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];
