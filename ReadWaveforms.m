(* ::Package:: *)

(* ::Chapter:: *)
(*Simulation Data Routines Package*)


(* ::Section:: *)
(*Begin ReadWaveforms Package*)


BeginPackage["KerrRingdownFitting`"]


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Section::Closed:: *)
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
