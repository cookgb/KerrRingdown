(* ::Package:: *)

(* ::Chapter:: *)
(*ReadWaveforms Package*)


(* ::Section:: *)
(*Begin ReadWaveforms Package*)


BeginPackage["ReadWaveforms`",{"DataRoutines`"}]


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Simulation Data Routines *)


ReadWaveforms::usage=
"ReadWaveforms[\!\(\*
StyleBox[\"dir\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"lm\",\nFontSlant->\"Italic\"]\),DataType\[Rule]\!\(\*
StyleBox[\"type\",\nFontSlant->\"Italic\"]\)] : Read in waveform data from a Numerical Relativity simulation.\n"<>
"\t \!\(\*
StyleBox[\"dir\",\nFontSlant->\"Italic\"]\) : String containing the full path to the directory containing simulation data\n"<>
"\t \!\(\*
StyleBox[\"lm\",\nFontSlant->\"Italic\"]\) : List of {l,m} pairs of spin-weighted spherical haromonic modes to import\n"<>
"\t \!\(\*
StyleBox[\"type\",\nFontSlant->\"Italic\"]\) : The type of data in file.  This option must be set.\n\n"<>
"Options:\n"<>
"\t DataType : Defaults to \!\(\*
StyleBox[\"None\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\".\",\nFontSlant->\"Italic\"]\)  Defined types: SXS\n"<>
"\t T0 : Defaults to 0.  Simulation time to set as T0 for ringdown fitting.\n"<>
"\t DataRange : Defaults to All.  Index range of time series data to read into memory.\n"<>
"\t RotateFrame : Defaults to False.  When it is set to True, a rotational transformation is done on the signal\n"<>
"\t to make the spin direction aligned with the QNM coordinates. When it is set to False, the orignal signal is used.\n\n"<>
"Additional options for each DataType. It can take the options from both SXSWaveform and SXSCCEWaveform. \n"<>
"\t SXS data:\n"<>
"\t\t SXSRNext : Defaults to 2.  Valid values 0, 2, 3, 4.\n"<>
"\t\t WaveformType: Defaults to Metric.  Valid values Psi4,Metric\n"<>
"\t\t FrameType : Defaults to Raw.  Valid values Raw,CoM."<>
"\t SXS_CCE data:\n"<>
"\t\t SXSRNext : Defaults to 2. CCE radii of a simulation or waveform extrapolation order to use.\n"<>
"\t\t WaveformType: Defaults to Metric.  Valid values Psi4, News, Metric\n"<>
"\t\t FrameType : Defaults to CoM.  Valid values CoM, Mem"<>
"\t Extrapolated : Defaults to False.\n"<>
"\t\t True : Extrapolated data is used. Read from {WaveformType}_Extrapolated_N{RNext}_CoM.h5\n"<>
"\t\t False : Extrapolated data is not used. CCE radii of a simulation is read in from RNext. Read from {WaveformType}_BondiCce_R{RNext}_CoM.h5\n"<>
"\t Superrest : Defaults to False.\n"<>
"\t\t True : Signal has been transformed to Super rest frame. Read from {WaveformType}..._CoM_Bondi.h5\n"<>
"\t\t False : Signal has not been tranformed to Super rest frame. Read from {Waveforminfo}..._CoM.h5\n"


TimeIndex::usage=
"TimeIndex[\!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\)] : Find the time index for a specific code time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) within the current series data stored in memeory."


PlotClm::usage=
"PlotClm[\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] : Returns a list of {t,\!\(\*SubscriptBox[\(C\), \(lm\)]\)} pairs suitable for plotting the full complex Numerical Relativity  "<>
"data for mode \!\(\*SubscriptBox[\(C\), \(lm\)]\)."


PlotReClm::usage=
"PlotReClm[\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] : Returns a list of {t,Re[\!\(\*SubscriptBox[\(C\), \(lm\)]\)]} pairs suitable for plotting the Real part of the Numerical Relativity  "<>
"data for mode \!\(\*SubscriptBox[\(C\), \(lm\)]\)."


PlotImClm::usage=
"PlotImClm[\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] : Returns a list of {t,Im[\!\(\*SubscriptBox[\(C\), \(lm\)]\)]} pairs suitable for plotting the Imaginary part of the Numerical Relativity  "<>
"data for mode \!\(\*SubscriptBox[\(C\), \(lm\)]\)."


PlotAbsClm::usage=
"PlotAbsClm[\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)] : Returns a list of {t,|\!\(\*SubscriptBox[\(C\), \(lm\)]\)|} pairs suitable for plotting the Magnitude of the Numerical Relativity  "<>
"data for mode \!\(\*SubscriptBox[\(C\), \(lm\)]\)."


PlotSumAbs2Clm::usage=
"PlotSumAbs2Clm[\!\(\*
StyleBox[\"lmlist\",\nFontSlant->\"Italic\"]\)] : Returns a list of {t,\!\(\*
StyleBox[\"\[CapitalSigma]\",\nFontSize->24]\)|\!\(\*SubscriptBox[\(C\), \(lm\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)} pairs suitable for plotting the Sum of the Square-Magnitudes of the Numerical Relativity  "<>
"data for a list of modes {\!\(\*SubscriptBox[\(C\), \(lm\)]\)}.  \!\(\*
StyleBox[\"lmlist\",\nFontSlant->\"Italic\"]\) has the form {{\!\(\*SubscriptBox[\(l\), \(1\)]\),\!\(\*SubscriptBox[\(m\), \(1\)]\)},{\!\(\*SubscriptBox[\(l\), \(2\)]\),\!\(\*SubscriptBox[\(m\), \(2\)]\)},...}"


(* ::Subsection:: *)
(*Reserved Globals*)


Protect[T0,DataType,SXS,SXSCCE,SXSRNext,RotateFrame];


Begin["`Private`"]


Protect[KRFtime,KRFC];


(* ::Section:: *)
(*Simulation Data Routines*)


Options[ReadWaveforms] = {T0->0,DataRange->All,DataType->None,SXSRNext->2,WaveformType->Metric,FrameType->Raw,RotateFrame->False,Extrapolated->False,Superrest->False};
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
						"BHPT",BHPTWaveform[fileordir,l,m,FilterRules[{opts},Options[BHPTWaveform]]],
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


TimeIndex[t_]:=If[t>=KRFtime[[-1]],Length[KRFtime],If[t<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=t},1][[1,1]]]]


PlotClm[l_,m_]:=Transpose[{KRFtime,KRFC[l,m]}]


PlotReClm[l_,m_]:=Transpose[{KRFtime,Re[KRFC[l,m]]}]


PlotImClm[l_,m_]:=Transpose[{KRFtime,Im[KRFC[l,m]]}]


PlotAbsClm[l_,m_]:=Transpose[{KRFtime,Abs[KRFC[l,m]]}]


PlotSumAbs2Clm[lm_List]:=Transpose[{KRFtime,Total[Function[x,Abs[KRFC[x[[1]],x[[2]]]]^2]/@lm]}]


(* ::Section::Closed:: *)
(*End of ReadWaveforms Package*)


End[] (* `Private` *)


EndPackage[]
