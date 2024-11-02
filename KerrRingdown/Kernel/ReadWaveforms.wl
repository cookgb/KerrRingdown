(* ::Package:: *)

(* ::Chapter:: *)
(*ReadWaveforms Package*)


(* ::Section::Closed:: *)
(*Begin ReadWaveforms Package*)


BeginPackage["KerrRingdown`"]


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Simulation Data Routines *)


ReadWaveforms::usage=
"ReadWaveforms[directory,{{\!\(\*SubscriptBox[\(l\), \(1\)]\), \!\(\*SubscriptBox[\(m\), \(1\)]\)},{\!\(\*SubscriptBox[\(l\), \(2\)]\), \!\(\*SubscriptBox[\(m\), \(2\)]\)},...},DataType\[RightArrow]Type] "<>
"Read in gravitational wave waveforms from a Numerical Relativity simulation file path "<>
"specified by directory.  DataType is an option that must be set."


TimeIndex::usage=
"TimeIndex[t] "<>
"Find the index for the currently stored time-series data that has a "<>
"time closest to the specified value of t."


PlotClm::usage=
"PlotClm[l,m] "<>
"Returns a list of {t,\!\(\*SubscriptBox[\(C\), \(lm\)]\)} pairs suitable for plotting the full "<>
"complex time-series data Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t)]."


PlotReClm::usage=
"PlotReClm[l,m] "<>
"Returns a list of {t,Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)]} pairs suitable for plotting the real part "<>
"of the full complex time-series data Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t)]."


PlotImClm::usage=
"PlotImClm[l,m] "<>
"Returns a list of {t,Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)]} pairs suitable for plotting the imaginary part "<>
"of the full complex time-series data Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t)]."


PlotAbsClm::usage=
"PlotAbsClm[l,m] "<>
"Returns a list of {t,Abs[\!\(\*SubscriptBox[\(C\), \(lm\)]\)]} pairs suitable for plotting the magnitude "<>
"of the full complex time-series data |\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t)|."


PlotSumAbs2Clm::usage=
"PlotSumAbs2Clm[{{\!\(\*SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},{\!\(\*SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)},...}] "<>
"Returns a list of {t, \!\(\*UnderscriptBox[\(\[Sum]\), \(lm\)]\)|\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)} pairs suitable for plotting the sum of the "<>
"square-magnitudes of the full complex time-series data \!\(\*UnderscriptBox[\(\[Sum]\), \(lm\)]\)|\!\(\*SubscriptBox[\(C\), \(lm\)]\)(t)\!\(\*SuperscriptBox[\(|\), \(2\)]\)."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


(*Protect[T0,DataType,SXS,SXSCCE,SXSRNext,RotateFrame];*)


Begin["`Private`"]


Protect[KRFtime,KRFC];


(* ::Section::Closed:: *)
(*Simulation Data Routines*)


Options[ReadWaveforms] = Union[{T0->0,DataType->None,SXSRNext->2,RotateFrame->False},Options[SXSWaveform],Options[SXSCCEWaveform]];
ReadWaveforms[fileordir_String,lm_List,opts:OptionsPattern[]]:= 
Module[{i,l,m,t,mode,lmodes,mmodes,t0=OptionValue[T0],
		range=OptionValue[DataRange],dattype=SymbolName[OptionValue[DataType]],
		sxsRnext=OptionValue[SXSRNext],tmpC,rot=False,tset=False},
	If[Head[OptionValue[RotateFrame]]==List,rot=True;{\[Theta],\[Phi]}=OptionValue[RotateFrame]];
	lmodes=DeleteDuplicates[#[[1]]&/@lm];
	Unprotect[KRFtime,KRFC];
	Clear[KRFtime,KRFC];
	Do[
		mmodes=Range[-l,l];
		Clear[tmpC];
		Do[
			If[rot || MemberQ[lm,{l,m}],
				{t,mode}=Switch[dattype,
						"SXS",SXSWaveform[fileordir,sxsRnext,l,m,FilterRules[{opts},Options[SXSWaveform]]],
						"SXSCCE",SXSCCEWaveform[fileordir,sxsRnext,l,m,FilterRules[{opts},Options[SXSCCEWaveform]]],
						_,Abort[]];
				If[!tset,KRFtime=t-t0;tset=True];
				If[!rot,KRFC[l,m]=mode;If[$MinPrecision>0,KRFC[l,m]=SetPrecision[KRFC[l,m],$MinPrecision]],tmpC[m]=mode];
			];
		,{m,mmodes}];
		If[rot,
			Do[
				If[MemberQ[lm,{l,m}],
					KRFC[l,m]=Sum[tmpC[mp]Conjugate[WignerD[{l,-mp,-m},\[Phi],\[Theta],0]],{mp,mmodes}];
					If[$MinPrecision>0,KRFC[l,m]=SetPrecision[KRFC[l,m],$MinPrecision]]
				];
			,{m,mmodes}];
		];
	,{l,lmodes}];
   Protect[KRFtime,KRFC];
]


(*Version 1: TimeIndex will find the next nearest time in KRFtime*)
(*TimeIndex[t_]:=If[t>=KRFtime[[-1]],Length[KRFtime],If[t<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=t},1][[1,1]]]]*)
(*Version 2: TimeIndex will find the nearest time in KRFtime*)
TimeIndex[t_]:=If[t>=KRFtime[[-1]],Length[KRFtime],If[t<=KRFtime[[1]],1,SequencePosition[KRFtime,Nearest[KRFtime,t,1],1][[1,1]]]]


PlotClm[l_,m_]:=Transpose[{KRFtime,KRFC[l,m]}]


PlotReClm[l_,m_]:=Transpose[{KRFtime,Re[KRFC[l,m]]}]


PlotImClm[l_,m_]:=Transpose[{KRFtime,Im[KRFC[l,m]]}]


PlotAbsClm[l_,m_]:=Transpose[{KRFtime,Abs[KRFC[l,m]]}]


PlotSumAbs2Clm[lm_List]:=Transpose[{KRFtime,Total[Function[x,Abs[KRFC[x[[1]],x[[2]]]]^2]/@lm]}]


(* ::Section::Closed:: *)
(*End of ReadWaveforms Package*)


End[] (* `Private` *)


EndPackage[]
