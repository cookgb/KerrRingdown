(* ::Package:: *)

(* ::Chapter:: *)
(*Simulation Data Routines Package*)


(* ::Section:: *)
(*Begin ReadWaveforms Package*)


BeginPackage["ReadWaveforms`"]


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
