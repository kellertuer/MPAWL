(* ::Package:: *)

(* ::Subsubtitle:: *)
(*The*)


(* ::Title:: *)
(*Multivariate periodic anisotropic Wavelet Library*)


(* ::Subtitle:: *)
(*Main File*)


(* ::Subsubtitle:: *)
(*Algorithms to perform the multivariate periodic Wavelet Analysis for anistropic Wavelets*)


(* ::Text:: *)
(*This Library can be included using Needs["MPAWL`"] and provides Functions for the multivariate periodic Analysis with anisotropic wavelets.*)


(* ::Program:: *)
(*Author: 		Ronny Bergmann*)
(*Created: 		13.11.2012*)
(*Last Changed: 	15.11.2012*)


BeginPackage["MPAWL`",
{
(*External dependencies*)
"SmithFormV6`" (* provided in this package, written by
Adriano Pascoletti, see
http://library.wolfram.com/infocenter/MathSource/7081/
*)
,
(*'internal' dependencies, i.e. sub packages *)
"MPAWL`Basics`",
"MPAWL`BoxSplines`",
"MPAWL`Pattern`",
"MPAWL`genSet`",
"MPAWL`TISpace`",
"MPAWL`Transforms`",
"MPAWL`Sampling`",
"MPAWL`deLaValleePoussin`",
"MPAWL`Visualization`"
}];


(* ::Section:: *)
(*Global Function Declaration*)


(* ::Section:: *)
(*Begin of private Area*)


Begin["`Private`"];


(* ::Subsection:: *)
(*End of private Function and Package Area*)


End[ ];


EndPackage[ ];