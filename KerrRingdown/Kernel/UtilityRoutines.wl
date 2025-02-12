(* ::Package:: *)

(* ::Chapter:: *)
(*Utility Routine Package*)


(* ::Section::Closed:: *)
(*Begin Utility Routine Package*)


BeginPackage["KerrRingdown`"]


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Utility Routines*)


QNModes::usage=
"QNModes[l,n] "<>
"Create a list of all valid QNM triplets for a single "<>
"value or a range of values of l, and for a single "<>
"value or range of overtones n.\n"<>
"QNModes[l,n] "<>
"Same as QNModes[l,n] except that the list is restricted "<>
"to the specified value or range of values for m."


SimulationModes::usage=
"SimulationModes[l] "<>
"Create a list of all possible  simulation harmonic modes for "<>
"a single value or a range of values of l.\n"<>
"SimulationModes[l,m] "<>
"Same as SimulationModes[l] except that the list is restricted "<>
"to the specified value or range of values for m."


QNMpIndex::Abort="Invalid QNM input. The input QNM should be a triplet {l, m, n}. "
QNMpIndex::usage=
"QNMpIndex[QNModesp,QNModesm,{l,m,n}] "<>
"Return the position index of the specified mode \!\(\*SubsuperscriptBox[\(C\), \(lmn\), \(+\)]\) in the "<>
"combined lists of QNMs QNModesp and QNModesm."


QNMmIndex::Abort="Invalid QNM input. The input QNM should be a triplet {l, m, n}. "
QNMmIndex::usage=
"QNMmIndex[QNModesp,QNModesm,{l,m,n}] "<>
"Return the position index of the specified mode \!\(\*SubsuperscriptBox[\(C\), \(lmn\), \(-\)]\) in the "<>
"combined lists of QNMs QNModesp and QNModesm."


SphericalHarmonicModes::Abort="Invalid QNM : `1`";
SphericalHarmonicModes::usage=
"SphericalHarmonicModes[{l,m,n},smodes] "<>
"Select the subset of signal modes in the list "<>
"\!\(\*StyleBox[\"smodes\", \"TI\"]\) that can overlap "<>
"with the QNM specified by {l,m,n}"


SphericalHarmonicModesNL::Abort="Invalid Quadratic Modes: `1`";
SphericalHarmonicModesNL::usage=""


SpheroidalHarmonicModes::Abort="Invalid simulation mode : `1`";
SpheroidalHarmonicModes::usage=
"SpheroidalHarmonicModes[{l,m},qnms] "<>
"Select the subset of QNMs in the list "<>
"\!\(\*StyleBox[\"qnms\", \"TI\"]\) that can overlap "<>
"with the signal mode specified by {l,m}."


SpheroidalHarmonicModesNL::Abort="Invalid simulation mode : `1`";
SpheroidalHarmonicModesNL::usage=""


FitMode::usage=
"FitMode[fit,{\!\(\*SubscriptBox[\(l\), \(s\)]\),\!\(\*SubscriptBox[\(m\), \(s\)]\)},\!\(\*SubscriptBox[\(t\), \(i\)]\)] "<>
"This function uses the output \!\(\*StyleBox[\"fit\", \"TI\"]\) from OverlapFit to "<>
"reconstruct the GW waveform corresponding to signal "<>
"mode \!\(\*FormBox[\(\*SubscriptBox[\(C\), \(\*SubscriptBox[\(l\), \(s\)] \*SubscriptBox[\(m\), \(s\)]\)](t)\),
TraditionalForm]\), based on the QNM expansion coefficients "<>
"fit at time \!\(\*SubscriptBox[\(t\), \(i\)]\)."


SimulateWaveforms::usage=
"SimulateWaveforms[t,SimModes,QNModesp,QNModesm,BHproperties] "<>
"Create a set of signal modes for a simulated ring-down "<>
"waveform based on a combination of specified QNMs."


QNMData::usage=
"QNMData[{l,m,n},l'] "<>
"Returns a list {\!\(\*SubscriptBox[\(\[Omega]\), \(lmn\)]\),\!\(\*SubscriptBox[\(A\), \(l' lmn\)]\)} containing the complex QNM "<>
"frequence and Spheroidal Harmonic expansion coefficient "<>
"for the current spin \[Chi] of the remnant black hole."


ViewRemnantParameterSpace\[Delta]\[Chi]::usage=
"ViewRemnantParameterSpace\[Delta]\[Chi][rps,\!\(\*SubscriptBox[\(t\), \(i\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \(i\)]\)] "<>
"produces a density plot and determines the parameters of the "<>
"minimum mismatch found in a 2-dimensional slice of the remnant "<>
"parameter space produced by RemnantParameterSearch."


ViewRemnantParameterSpace::usage=
"ViewRemnantParameterSpace[rps] "<>
"Interactively view  the results from RemnantParameterSearch."


Add\[Delta]\[Chi]Lines::usage=
"Add\[Delta]\[Chi]Lines[\[Delta]\[Chi]plot] "<>
"adds lines to a density plot to help visualizing the location "<>
"of the minimum mismatch."


Add\[Delta]\[Chi]LinesandLabel::usage=
"Add\[Delta]\[Chi]LinesandLabel[\[Delta]\[Chi]plot,label] "<>
"adds lines to a density plot to help visualizing the location "<>
"of the minimum mismatch, and include a label to display "<>
"additional information."


OverlapSequenceAmplitudes::badmodel="Model information in fitinfo has wrong length";
OverlapSequenceAmplitudes::usage=
"OverlapSequenceAmplitudes[fit] "<>
"Returns the QNM expansion coefficients and their standard errors for "<>
"each mode from an OverlapFit result "<>
"\!\(\*StyleBox[\"fit\", \"TI\"]\). The results are presented as amplitudes and phases."


OverlapSequenceCoefPlus::Abort="QNM input does not exist in OverlapFit result \!\(\*
StyleBox[\"fit\",\nFontSlant->\"Italic\"]\). "
OverlapSequenceCoefPlus::usage=
"OverlapSequenceCoefPlus[fit,{l,m,n}] "<>
"Returns the time sequence of the complex QNM expansion coefficient "<>
"\!\(\*SubsuperscriptBox[\(C\), \(lmn\), \(+\)]\)(t) from an OverlapFit result fit"


OverlapSequenceCoefMinus::Abort="QNM input does not exist in OverlapFit result \!\(\*
StyleBox[\"fit\",\nFontSlant->\"Italic\"]\). "
OverlapSequenceCoefMinus::usage=
"OverlapSequenceCoefMinus[fit,{l,m,n}] "<>
"Returns the time sequence of the complex QNM expansion coefficient "<>
"\!\(\*SubsuperscriptBox[\(C\), \(lmn\), \(-\)]\)(t) from an OverlapFit result fit"


OverlapSequenceCoefNL::Abort="Quadratic mode input does not exist in OverlapFit result \!\(\*
StyleBox[\"fit\",\nFontSlant->\"Italic\"]\). "
OverlapSequenceCoefMinus::usage=""


MOSAmp::usage=
"MOSAmp[ampdata,index] "<>
"returns a list appropriate for plotting QNM amplitudes as a "<>
"function of time for the QNM specified by index.  The data "<>
"is taken from ampdat as returned by "<>
"MaxOverlapSequenceAmplitudes."


MOSPhase::usage=
"MOSPhase[ampdata,index] "<>
"returns a list appropriate for plotting scaled QNM phases "<>
"(\[Phi]/\[Pi]) as a function of time for the QNM specified by index.  "<>
"The data is taken from ampdat as returned by "<>
"MaxOverlapSequenceAmplitudes."


MergeRPSLists::parametermismatch="Two RPS lists being merged have different parameter ranges or use different sets of modes.";
MergeRPSLists::usage=
"MergeRPSLists[\!\(\*SubscriptBox[\(rps\), \(1\)]\),\!\(\*SubscriptBox[\(rps\), \(2\)]\),...] "<>
"combines two or more ring-down parameter search lists "<>
"\!\(\*SubscriptBox[\(rps\), \(i\)]\) as returned by RemnantParameterSearch. "<>
"The returned list is sorted in time and duplicate time "<>
"entries are removed."


MergeMaxOverlapSequences::usage=
"MergeMaxOverlapSequences[\!\(\*SubscriptBox[\(mos\), \(1\)]\),\!\(\*SubscriptBox[\(mos\), \(2\)]\),...] "<>
"combines two or more maximum overlap sequence lists "<>
"\!\(\*SubscriptBox[\(mos\), \(i\)]\) as returned by RemnantParameterSpaceMaxOverlap.  "<>
"The returned list is sorted in time and duplicate time "<>
"entries are removed."


FindNonlinearIndex::usage=""


SetGreedyModes::usage=""


FitInfoStruct::usage=""


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*Utility Routines*)


QNModes[la_Integer|la_List,na_Integer|na_List]:=
Module[{s=-2,li,ni,l,m,n,llist,nlist,modes={}},
   llist=Cases[DeleteDuplicates[Sort[If[Head[la]==Integer,{la},Null[],la]]],x_/;x>=Abs[s]->x];
   nlist=Cases[DeleteDuplicates[Sort[If[Head[na]==Integer,{na},Null[],na]]],x_/;x>=0->x];
   For[li=1,li<=Length[llist],++li,l=llist[[li]];
      For[ni=1,ni<=Length[nlist],++ni,n=nlist[[ni]];
         For[m=l,m>=-l,--m,AppendTo[modes,{l,m,n}]
         ];
      ];
   ];
   modes
]


QNModes[la_Integer|la_List,ma_Integer|ma_List,na_Integer|na_List]:=
Module[{modes},
	modes=QNModes[la,na];
	Switch[Head[ma],
		List,Cases[modes,{_,m_,_}/;MemberQ[ma,m]],
		Integer,Cases[modes,{_,ma,_}]
	]
]


SimulationModes[la_Integer|la_List]:=
Module[{s=-2,li,l,m,llist,modes={}},
   llist=Cases[DeleteDuplicates[Sort[If[Head[la]==Integer,{la},Null[],la]]],x_/;x>=Abs[s]->x];
   For[li=1,li<=Length[llist],++li,l=llist[[li]];
      For[m=l,m>=-l,--m,AppendTo[modes,{l,m}]
      ];
   ];
   modes
]


SimulationModes[la_Integer|la_List,ma_Integer|ma_List]:=
Module[{modes},
	modes=SimulationModes[la];
	Switch[Head[ma],
		List,Cases[modes,{_,m_}/;MemberQ[ma,m]],
		Integer,Cases[modes,{_,ma}]
	]
]


QNMpIndex[QNModesp_List,QNModesm_List,qnm_List]:=
Module[{p},
	If[Dimensions[qnm]!={3},Message[QNMpIndex::Abort];Abort[]];
	p=Position[QNModesp,qnm];
	If[Length[p]==0,0,p[[1,1]]]
]


QNMmIndex[QNModesp_List,QNModesm_List,qnm_List]:=
Module[{p},
	If[Dimensions[qnm]!={3},Message[QNMmIndex::Abort];Abort[]];
	p=Position[QNModesm,qnm];
	If[Length[p]==0,0,Length[QNModesp]+p[[1,1]]]
]


FindNonlinearIndex[QNModesp_List,QNModesm_List,NLlist_List,NLfind_List]:=Module[{p},
p=Position[NLlist,NLfind];
If[Length[p]==0,0,Length[QNModesp]+Length[QNModesm]+p[[1,1]]]]


SphericalHarmonicModes[qnm_List,sim_List]:=
Module[{s=-2},
   If[Length[qnm]!=3||qnm[[1]]<Abs[s]||Abs[qnm[[2]]]>qnm[[1]]||qnm[[3]]<0,
      Message[SphericalHarmonicModes::Abort,qnm];Abort[]];
   DeleteDuplicates[Cases[sim,x_/;x[[1]]>=Max[Abs[s],Abs[qnm[[2]]]]->x]]
]


(*This function works the same as the SphericalHarmonicModes but for quadratic modes. *)
(*qqnm is a list in form of {l,m,\[Omega]}*)
SphericalHarmonicModesNL[qqnm_List,sim_List]:=
Module[{s=-2},
   If[Length[qqnm]!=3,
      Message[SphericalHarmonicModesNL::Abort,qqnm];Abort[]];
   DeleteDuplicates[Cases[sim,x_/;x[[1]]>=Max[Abs[s],Abs[qqnm[[2]]]]->x]]
]


SpheroidalHarmonicModes[sim_List,qnm_List]:=
Module[{s=-2},
   If[Length[sim]!=2||sim[[1]]<Abs[s]||Abs[sim[[2]]]>sim[[1]],
      Message[SpheroidalHarmonicModes::Abort,sim];Abort[]];
   DeleteDuplicates[Cases[Cases[qnm,x_/;x[[1]]>=Max[Abs[s],Abs[sim[[2]]]]->x],x_/;Abs[x[[2]]]<=sim[[1]]->x]]
]


(*This function works the same as the SpheroidalHarmonicModes but for quadratic modes. *)
(*qqnm is a list in form of {l,m,\[Omega]}*)
SpheroidalHarmonicModesNL[sim_List,qqnm_List]:=
Module[{s=-2},
	If[Length[sim]!=2||sim[[1]]<Abs[s]||Abs[sim[[2]]]>sim[[1]],
      Message[SpheroidalHarmonicModesNL::Abort,sim];Abort[]];
    DeleteDuplicates[Cases[Cases[qqnm,x_/;x[[1]]>=Max[Abs[s],Abs[sim[[2]]]]->x],x_/;Abs[x[[2]]]<=sim[[1]]->x]] 
]


Options[FitMode]= Union[Options[SetModeData],{OmitModes->{{},{},{}}}]
FitMode[fit_List,lm_List,fittime_?NumberQ,opts:OptionsPattern[]]:=
Module[{t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,err2,count,qnmp,qnmm,l,m,lp,mp,np,ind,
        j,mode=Table[0,Length[KRFtime]],omitp,omitm,
        NLmodesExist,omitNL,
        fixedGreedyIndex={},FixedModesG={},QQNModes,
		qnmpSet,qnmmSet,QNModespUnfixed,QNModesmUnfixed,massratioFixed,QNModespFixed,
		QNModesmFixed,QQNModesFixed={},QQNModesUnfixed,QQNModes\[Omega]Unfixed,QQNModes\[Omega]Fixed,
		\[Omega]NL,k,fixedAmp,qqnm,qqnmIndex},
	{omitp,omitm,omitNL}=OptionValue[OmitModes];
	{l,m}=lm;
	(* obtain the information of QNMp, QNMm, Quadratic QNM, and fixed modes. 
	Added for making this function compatible with qudratic QNMs and greedy algorithm.
	*) 	
	{{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{QNModesp,QNModesm,QQNModes},{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}} = FitInfoStruct[fit];
	If[QQNModes=={},NLmodesExist=False,NLmodesExist=True,NLmodesExist=True];	
	SetModeData[a,qnmpSet,qnmmSet,Evaluate@FilterRules[{opts},Options@SetModeData]];
	If[Length[FixedModesG]!=0,
		massratioFixed = fixedGreedyIndex[[-1]];
		QNModespFixed=QNModesp[[#]]&/@fixedGreedyIndex[[1]];
		QNModesmFixed=QNModesm[[#]]&/@fixedGreedyIndex[[2]];
		QQNModesFixed=QQNModes[[#]]&/@fixedGreedyIndex[[3]];
	];
	If[NLmodesExist,(*calculate the frequency of the quadratic QNMs*)
		QQNModes\[Omega]Unfixed=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@QQNModesUnfixed];
		QQNModes\[Omega]Fixed=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@QQNModesFixed];
	];
	
	ind=If[fittime>=t[[-1]],Length[t],If[fittime<=t[[1]],1,SequencePosition[t,{x_/;x>=fittime},1][[1,1]]]];
	
	(*contribution to waveform from unfixed QNM plus*)
	qnmp=SpheroidalHarmonicModes[{l,m},QNModespUnfixed];
	For[j=1,j<=Length[qnmp],++j,
		If[MemberQ[omitp,qnmp[[j]]],Continue[]];
		{lp,mp,np}=qnmp[[j]];
		mode+=If[UseSpheroidalExpansion,
				WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
				(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
				KroneckerDelta[l,lp]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
		];
	];
	(*contribution to waveform from unfixed QNM minus*)
	qnmm=SpheroidalHarmonicModes[{l,m},QNModesmUnfixed];
	For[j=1,j<=Length[qnmm],++j,
		If[MemberQ[omitm,qnmm[[j]]],Continue[]];
		{lp,mp,np}=qnmm[[j]];
		mode+=If[UseSpheroidalExpansion,
			(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
			KroneckerDelta[l,lp]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
		];
	];
	If[NLmodesExist, (*if there are quadratic QNMs*)
		qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Unfixed];
		qqnmIndex=Flatten[Position[QQNModes\[Omega]Unfixed,#]&/@qqnm];
		For[j=1,j<=Length[qqnm],j++,
			If[MemberQ[omitNL,QQNModesUnfixed[[qqnmIndex[[j]]]]],Continue[]];
			{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
			mode+=WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*amp[[ind,FindNonlinearIndex[QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,QQNModesUnfixed[[qqnmIndex[[j]]]]]]]*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* KRFtime/massratio]	
		];
	];
	If[Length[FixedModesG]!=0,(*if there are fixed QNMs in the fit results*)
			(*contribution to waveform from fixed QNM plus and minus*)
			qnmp=SpheroidalHarmonicModes[{l,m},QNModespFixed];
			For[j=1,j<=Length[qnmp],++j,
				If[MemberQ[omitp,qnmp[[j]]],Continue[]];
				{lp,mp,np}=qnmp[[j]];
				For[k=1,k<=Length[FixedModesG],k++,
					If[FixedModesG[[k,1]]=={lp,mp,np,1},fixedAmp=FixedModesG[[k,2]];Break[]]
				];
				mode+=If[UseSpheroidalExpansion,
					WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratioFixed],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]*fixedAmp*Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratioFixed]
				];
			];
			qnmm=SpheroidalHarmonicModes[{l,m},QNModesmFixed];
			For[j=1,j<=Length[qnmm],++j,
				If[MemberQ[omitm,qnmm[[j]]],Continue[]];
				{lp,mp,np}=qnmm[[j]];
				For[k=1,k<=Length[FixedModesG],k++,
					If[FixedModesG[[k,1]]=={lp,mp,np,-1},fixedAmp=FixedModesG[[k,2]];Break[]]
				];
				mode+=If[UseSpheroidalExpansion,
						(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratioFixed],
						(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
						KroneckerDelta[l,lp]*fixedAmp*Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratioFixed]
					];
			];
			If[NLmodesExist&&QQNModes\[Omega]Fixed!={}, (*if there are quadratic QNMs*)
				qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Fixed];
				qqnmIndex=Flatten[Position[QQNModes\[Omega]Fixed,#]&/@qqnm];
				For[j=1,j<=Length[qqnm],j++,
					If[MemberQ[omitNL,QQNModesFixed[[qqnmIndex[[j]]]]],Continue[]];
					{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
					For[k=1,k<=Length[FixedModesG],k++,
						If[FixedModesG[[k,1]]==QQNModesFixed[[qqnmIndex[[j]]]],fixedAmp=FixedModesG[[k,2]];Break[]]
					];
					mode+=WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* KRFtime/massratioFixed]	
				];
			];
	];
	{KRFtime,mode}
]


Options[SimulateWaveforms]=Union[{Width->5,T0->0},Options[SetModeData]];
SimulateWaveforms[t_List,SimModes_List,QNModesp_List,QNModesm_List,BHparameters_List,opts:OptionsPattern[]]:=
Module[{massratio,a,\[Theta],\[Phi],qnmp,qnmm,amp,i,l,m,fitinfo,envelope,width=OptionValue[Width],t0=OptionValue[T0]},
   {massratio,a,\[Theta],\[Phi]}=BHparameters;
   qnmp=#[[1]]&/@QNModesp;qnmm=#[[1]]&/@QNModesm;
   amp={Join[{1,I} . FromPolarCoordinates[#[[2]]]&/@QNModesp,{1,I} . FromPolarCoordinates[#[[2]]]&/@QNModesm]};
   SetModeData[a,qnmp,qnmm,Evaluate@FilterRules[{opts},Options@SetModeData]];
   fitinfo={t,Null[],amp,{BHparameters,SimModes,qnmp,qnmm},Null[],Null[]};
   envelope=(ArcTan[10(t-t0)/width]+\[Pi]/2)/\[Pi];
   Unprotect[KRFtime,KRFC];
   Clear[KRFtime,KRFC];
   KRFtime=t;
   For[i=1,i<=Length[SimModes],++i,
      {l,m}=SimModes[[i]];
      KRFC[l,m]= envelope(FitMode[fitinfo,{l,m},t[[1]]][[2]]);
   ];
   Protect[KRFtime,KRFC];
]


QNMData[qnm_List,lp_]:={KRF\[Omega][qnm[[1]],qnm[[2]],qnm[[3]]],KRFYS[lp,qnm[[1]],qnm[[2]],qnm[[3]]]}


Options[ViewRemnantParameterSpace\[Delta]\[Chi]]={Log10MisMatchRange->Automatic};
ViewRemnantParameterSpace\[Delta]\[Chi][rps_,ti_,\[Theta]i_,OptionsPattern[]]:=
Module[{time,mismatch,amp,BHparams,SimModes,QNModesp,QNModesm,err,count,data,min,
		range=OptionValue[Log10MisMatchRange],colorfunc,rng,NLTemp,GreedyTemp},
	Switch[Length[rps[[4]]],
		4,
		{time,mismatch,amp,{BHparams,SimModes,QNModesp,QNModesm},err,count}=rps,
		5,
		{time,mismatch,amp,{BHparams,SimModes,QNModesp,QNModesm,NLTemp},err,count}=rps,
		6,
		{time,mismatch,amp,{BHparams,SimModes,QNModesp,QNModesm,NLTemp,GreedyTemp},err,count}=rps
	];	
	(*data = Transpose[{Flatten[Take[BHparams,All,All,{\[Theta]i},1,{2}]],
						Flatten[Take[BHparams,All,All,{\[Theta]i},1,{1}]],
						Log10[Flatten[Take[mismatch,{ti},All,All,{\[Theta]i}]]]}];*)
	data = Transpose[{Flatten[Take[BHparams,All,All,{\[Theta]i},{2}]],
						Flatten[Take[BHparams,All,All,{\[Theta]i},{1}]],
						Log10[Flatten[Take[mismatch,{ti},All,All,{\[Theta]i}]]]}];
	min=Flatten[MinimalBy[data,#[[3]]&]];
	colorfunc=Blend[Transpose[{{0,0.05,0.1,0.2,0.4,0.6,0.8,1.0},{White,Lighter[Lighter[Orange]],Lighter[Orange],Orange,Red,Darker[Red],Darker[Darker[Red]],Black}}],#1]&;
	If[Head[range]==List,
		colorfunc=Blend[Transpose[{{0,0.05,0.1,0.2,0.4,0.6,0.8,1.0},{White,Lighter[Lighter[Orange]],Lighter[Orange],Orange,Red,Darker[Red],Darker[Darker[Red]],Black}}],Rescale[#1,rng]]&/.rng->range;
	];
	{Labeled[
		ListDensityPlot[data,
						PlotRange->{Full,Full,range},
						PlotLegends->If[Head[range]==List,BarLegend[{Automatic,range}],Automatic,Automatic],
						ColorFunctionScaling->If[Head[range]==List,False,True,True],
						ColorFunction->colorfunc,
						LabelStyle->Directive[Black,16],PlotRangePadding->None,ImageSize->400],
			{Style["\!\(\*
StyleBox[SubscriptBox[
StyleBox[\"M\",\nFontWeight->\"Plain\"], \"f\"],\nFontWeight->\"Plain\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"[\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"M\",\nFontWeight->\"Plain\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"]\",\nFontWeight->\"Plain\"]\)",20],Style["\!\(\*
StyleBox[SubscriptBox[\"\[Chi]\", 
StyleBox[\"f\",\nFontSlant->\"Italic\"]],\nFontWeight->\"Plain\"]\)",20],Style["\!\(\*
StyleBox[SubscriptBox[\"log\", \"10\"],\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"\[ScriptCapitalM]\",\nFontWeight->\"Plain\"]\)",20]},{Left,Bottom,Right},RotateLabel->True],min}
]


ViewRemnantParameterSpace[rps_]:=
DynamicModule[{ti=1,\[Theta]i=1,range=Automatic,save={{0,0,0},Null[]}},
	Panel[
		(*Column[{
			Row[{Slider[Dynamic[ti],{1,Length[rps[[1]]],1},ImageSize->350],
				Dynamic[TextCell[Style["  t = "<>ToString[rps[[1,ti]]],14]]]}],
			Row[{Slider[Dynamic[\[Theta]i],{1,Dimensions[rps[[4,1]]][[3]],1},ImageSize->350],
				Dynamic[TextCell[Style["  \[Theta] = "<>ToString[Take[rps[[4,1]],All,All,{\[Theta]i}][[1,1,1,1,3]],FormatType->InputForm],14]]]}],
			Row[{Dynamic[TextCell[Style["Minimum : "<>ToString[save[[2]]],14]]]}],
			Row[{Dynamic[ControlActive[Show[save[[1]]],save=ViewRemnantParameterSpace\[Delta]\[Chi][rps,ti,\[Theta]i];save[[1]]]]}]
			}]*)
		Column[{
			Row[{Slider[Dynamic[ti],{1,Length[rps[[1]]],1},ImageSize->350],
				Dynamic[TextCell[Style["  t = "<>ToString[rps[[1,ti]]],14]]]}],
			Row[{Slider[Dynamic[\[Theta]i],{1,Dimensions[rps[[4,1]]][[3]],1},ImageSize->350],
				Dynamic[TextCell[Style["  \[Theta] = "<>ToString[Take[rps[[4,1]],All,All,{\[Theta]i}][[1,1,1,3]],FormatType->InputForm],14]]]}],
			Row[{Dynamic[TextCell[Style["Minimum : "<>ToString[save[[2]]],14]]],
			TextCell[Style[" | DataRange ",14]],
			InputField[Dynamic[range],FieldSize->Small,BaseStyle->14]}],
			Row[{Dynamic[ControlActive[Show[save[[1]]],save=ViewRemnantParameterSpace\[Delta]\[Chi][rps,ti,\[Theta]i,Log10MisMatchRange->range];save[[1]]]]}]
			}]
		]
]


Add\[Delta]\[Chi]Lines[lp_]:=ReplacePart[lp,{1,1,1}->Show[lp[[1,1,1]],
				Graphics[{AbsoluteThickness[2],Dashed,White,Line[{{lp[[2,1]],0.85},{lp[[2,1]],lp[[2,2]]}}]}],
				Graphics[{AbsoluteThickness[2],Dashed,White,Line[{{0.5,lp[[2,2]]},{lp[[2,1]],lp[[2,2]]}}]}]]]


Add\[Delta]\[Chi]LinesandLabel[lp_,label_]:=ReplacePart[lp,{1,1,1}->Show[lp[[1,1,1]],
				Graphics[{AbsoluteThickness[2],Dashed,White,Line[{{lp[[2,1]],0.85},{lp[[2,1]],lp[[2,2]]}}]}],
				Graphics[{AbsoluteThickness[2],Dashed,White,Line[{{0.5,lp[[2,2]]},{lp[[2,1]],lp[[2,2]]}}]}],
				Graphics[label]]]


Options[OverlapSequenceAmplitudes]=Union[{},Options[SetModeData]];
OverlapSequenceAmplitudes[fit_List,opts:OptionsPattern[]]:=
Module[{fittime,t,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,
		ind,ind2,\[Rho],amp,err2,count,i,l,m,mode,qnmp,qnmm,j,lp,mp,np,
		cerr,err2i,amps={},err2s={},NLmodesExist,
		qnmTemp,fixedGreedyIndex={},FixedModesG={},QQNModes,
		qnmpSet,qnmmSet,QNModespUnfixed,QNModesmUnfixed,massratioFixed,QNModespFixed,
		QNModesmFixed,QQNModesFixed={},QQNModesUnfixed={},QQNModes\[Omega]Unfixed,QQNModes\[Omega]Fixed,
		\[Omega]NL,k,fixedAmp,qqnm,qqnmIndex},

	(* obtain the information of QNMp, QNMm, Quadratic QNM, and fixed modes. 
	Added for making this function compatible with qudratic QNMs and greedy algorithm.
	*)
	{{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{QNModesp,QNModesm,QQNModes},{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}} = FitInfoStruct[fit];
	If[QQNModes=={},NLmodesExist=False,NLmodesExist=True,NLmodesExist=True];
	SetModeData[a,qnmpSet,qnmmSet,Evaluate@FilterRules[{opts},Options@SetModeData]];
	If[Length[FixedModesG]!=0,
		massratioFixed = fixedGreedyIndex[[-1]];
		QNModespFixed=QNModesp[[#]]&/@fixedGreedyIndex[[1]];
		QNModesmFixed=QNModesm[[#]]&/@fixedGreedyIndex[[2]];
		QQNModesFixed=QQNModes[[#]]&/@fixedGreedyIndex[[3]];
	];
	If[NLmodesExist,
		QQNModes\[Omega]Unfixed=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@QQNModesUnfixed];
		QQNModes\[Omega]Fixed=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@QQNModesFixed];
	];
	For[ind=1,ind<=Length[t],++ind,
		fittime=t[[ind]];
		ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
		mode=Table[0,Length[SimModes],Length[KRFtime]];
		cerr=0;
		For[i=1,i<=Length[SimModes],++i,
			{l,m}=SimModes[[i]];
			mode=Table[0,Length[KRFtime]];
			(*contribution to waveform from unfixed QNM plus*)
			qnmp=SpheroidalHarmonicModes[{l,m},QNModespUnfixed];
			For[j=1,j<=Length[qnmp],++j,
				{lp,mp,np}=qnmp[[j]];
				mode+=If[UseSpheroidalExpansion,
					WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
				];
			];
			(*contribution to waveform from unfixed QNM minus*)
			qnmm=SpheroidalHarmonicModes[{l,m},QNModesmUnfixed];
			For[j=1,j<=Length[qnmm],++j,
				{lp,mp,np}=qnmm[[j]];
				mode+=If[UseSpheroidalExpansion,
					(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
				];
			];
			If[NLmodesExist, (*if there are quadratic QNMs*)
				qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Unfixed];
				qqnmIndex=Flatten[Position[QQNModes\[Omega]Unfixed,#]&/@qqnm];
				For[j=1,j<=Length[qqnm],j++,
					{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
					mode+=WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*amp[[ind,FindNonlinearIndex[QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,QQNModesUnfixed[[qqnmIndex[[j]]]]]]]*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* KRFtime/massratio]	
				];
			];
			If[Length[FixedModesG]!=0,(*if there are fixed QNMs in the fit results*)
					(*contribution to waveform from fixed QNM plus and minus*)
					qnmp=SpheroidalHarmonicModes[{l,m},QNModespFixed];
					For[j=1,j<=Length[qnmp],++j,
						{lp,mp,np}=qnmp[[j]];
						For[k=1,k<=Length[FixedModesG],k++,
							If[FixedModesG[[k,1]]=={lp,mp,np,1},fixedAmp=FixedModesG[[k,2]];Break[]]
						];
						mode+=If[UseSpheroidalExpansion,
						WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratioFixed],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
						KroneckerDelta[l,lp]*fixedAmp*Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratioFixed]
						];
					];
					qnmm=SpheroidalHarmonicModes[{l,m},QNModesmFixed];
					For[j=1,j<=Length[qnmm],++j,
						{lp,mp,np}=qnmm[[j]];
						For[k=1,k<=Length[FixedModesG],k++,
							If[FixedModesG[[k,1]]=={lp,mp,np,-1},fixedAmp=FixedModesG[[k,2]];Break[]]
						];
						mode+=If[UseSpheroidalExpansion,
						(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratioFixed],
						(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
						KroneckerDelta[l,lp]*fixedAmp*Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratioFixed]
						];
					];
					If[NLmodesExist&&QQNModes\[Omega]Fixed!={}, (*if there are quadratic QNMs*)
						qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Fixed];
						qqnmIndex=Flatten[Position[QQNModes\[Omega]Fixed,#]&/@qqnm];
						For[j=1,j<=Length[qqnm],j++,
							{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
							For[k=1,k<=Length[FixedModesG],k++,
								If[FixedModesG[[k,1]]==QQNModesFixed[[qqnmIndex[[j]]]],fixedAmp=FixedModesG[[k,2]];Break[]]
							];
							mode+=WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*fixedAmp*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* KRFtime/massratioFixed]	
						];
					];
			];
			mode-=KRFC[l,m];
        	mode=Take[Drop[mode,ind2-1],count[[ind]]];
        	cerr+=Conjugate[mode] . mode;
        ];
        err2i=Sqrt[err2[[ind]]*Abs[cerr]/(2(Length[SimModes]count[[ind]]-Length[QNModespUnfixed]-Length[QNModesmUnfixed]-Length[QQNModesUnfixed]))];
		AppendTo[amps,#/{1,\[Pi]}&/@AbsArg[amp[[ind]]]];
		AppendTo[err2s,Transpose[{err2i,If[Im[#]!=0,1,#/\[Pi]]&/@ArcSin[err2i/Abs[amp[[ind]]]]}]];
	];
	{t,Transpose[amps,2<->1],Transpose[err2s,2<->1]}
]


OverlapSequenceCoefPlus[fit_List,qnmp_List]:=Module[
{qnmpIndex},
	qnmpIndex = QNMpIndex[fit[[4,3]],fit[[4,4]],qnmp];
	If[qnmpIndex==0,Message[OverlapSequenceCoefPlus::Abort];Abort[]];
	Transpose[fit[[3]]][[qnmpIndex]]
]


OverlapSequenceCoefMinus[fit_List,qnmm_List]:=Module[
{qnmmIndex},
	qnmmIndex = QNMmIndex[fit[[4,3]],fit[[4,4]],qnmm];
	If[qnmmIndex==0,Message[OverlapSequenceCoefMinus::Abort];Abort[]];
	Transpose[fit[[3]]][[qnmmIndex]]
]


(*This is the quadratic mode's version of OverlapSequenceCoefPlus*)
OverlapSequenceCoefNL[fit_List,qqnm_List]:=Module[
{qqnmIndex},
	If[Length[fit[[4]]]<5||Length[Dimensions[fit[[4,5]]]]==1,Message[OverlapSequenceCoefNL::Abort];Abort[]];
	qqnmIndex = FindNonlinearIndex[fit[[4,3]],fit[[4,4]],fit[[4,5]],qqnm];
	If[qqnmIndex==0,Message[OverlapSequenceCoefNL::Abort];Abort[]];
	Transpose[fit[[3]]][[qqnmIndex]]
]


Options[MOSAmp]={Errors->False,Log10->False};
MOSAmp[ampdat_List,ampind_Integer,opts:OptionsPattern[]]:=
Module[{ploterrs=OptionValue[Errors],plotlog=OptionValue[Log10],t,amp,err},
	t = ampdat[[1]];
	amp = Flatten[Take[ampdat[[2]],{ampind},All,{1}]];
	If[ploterrs,
		err = Flatten[Take[ampdat[[3]],{ampind},All,{1}]];
		amp = Around[#[[1]],#[[2]]]&/@Transpose[{amp,err}];
	];
	If[plotlog,
		amp = Log10[#]&/@amp;
	];
	Transpose[{t,amp}]
]


Options[MOSPhase]={Errors->False};
MOSPhase[ampdat_List,ampind_Integer,opts:OptionsPattern[]]:=
Module[{ploterrs=OptionValue[Errors],t,amp,err},
	t = ampdat[[1]];
	amp = Flatten[Take[ampdat[[2]],{ampind},All,{2}]];
	If[ploterrs,
		err = Flatten[Take[ampdat[[3]],{ampind},All,{2}]];
		amp = Around[#[[1]],#[[2]]]&/@Transpose[{amp,err}];
	];
	Transpose[{t,amp}]
]


MergeRPSLists[rps1_List,rps2_List]:=Module[
{t,tordering,orderedB,dups,reorderB,i,mm,amp,bhp,err,cnt},
	If[rps1[[4]]!=rps2[[4]],Message[MergeRPSLists::parametermismatch];Abort[]];
	bhp=rps1[[4]];
	t=Join[rps1[[1]],rps2[[1]]];
	tordering=Ordering[t];
	dups=Flatten[Last@Reap[Sow[Drop[#,1]]&/@Select[GatherBy[Range@Length[t],t[[#]]&],Length[#]>1&]]];
	tordering=DeleteCases[tordering,n_/;MemberQ[dups,n]];
	orderedB=OrderedQ[tordering];
	reorderB= Not[orderedB] || Length[dups]>0;
	If[reorderB,t=t[[tordering]]];
	mm=If[reorderB,Join[rps1[[2]],rps2[[2]]][[tordering]],Join[rps1[[2]],rps2[[2]]]];
	amp=If[reorderB,Join[rps1[[3]],rps2[[3]]][[tordering]],Join[rps1[[3]],rps2[[3]]]];
	err=If[reorderB,Join[rps1[[5]],rps2[[5]]][[tordering]],Join[rps1[[5]],rps2[[5]]]];
	cnt=If[reorderB,Join[rps1[[6]],rps2[[6]]][[tordering]],Join[rps1[[6]],rps2[[6]]]];
	{t,mm,amp,bhp,err,cnt}
]
MergeRPSLists[n__,rps1_List,rps2_List]:=MergeRPSLists[n,MergeRPSLists[rps1,rps2]]


MergeMaxOverlapSequences[mol1_List,mol2_List]:=
	DeleteDuplicatesBy[SortBy[Join[mol1,mol2],First],First]
MergeMaxOverlapSequences[n__,mol1_List,mol2_List]:=
	MergeMaxOverlapSequences[n,MergeMaxOverlapSequences[mol1,mol2]]


SetGreedyModes[qnmpIn_List,qnmmIn_List,NLlist_,FixedModesG_List,param_]:=Module[
{qnmp=qnmpIn,qnmm=qnmmIn,FixedModes,FixedQNM={},FixedNL={},idxFixedP={},idxFixedM={},idxFixedNL={},FixedQNMp={},FixedQNMm={},i,FixedQnmParamP,FixedQnmParamM,
FixedNlParamP,FixedNlParamM,idxFixedNlP={},idxFixedNlM={},qnmpFixedNL={},qnmmFixedNL={},qnmpNL={},qnmmNL={},qnmpSet,qnmmSet,
FixedModesGPM={},FixedModesGP={},FixedModesGM={},FixedModesGNL={},FixedModesGsorted},
	(*
	This is a function called by the OverlapFit function to sort out the fixed QNMs from unknown QNMs
	Limitation: It can only deal with the case that the nonlinear mode pair have no overlap with the QNMs. 
	For example, we fit quadratic pair {{2,2,0,+},{2,0,0,+}} that contributes to spherical harmonics (2,2).
	When QNM {2,2,0,+} is also in the fitting set and to be fixed, it can cause problems. 
	To Do: Solve that limitation issue. 
	*)
If[Head[NLlist]==List,If[#[[4]]==1,AppendTo[qnmpNL,#[[1;;3]]],AppendTo[qnmmNL,#[[1;;3]]]]&/@DeleteDuplicates[Flatten[NLlist,1]]];
(*Separate the fixed QNMs and the fixed nonlinear modes*)
	FixedModes=#[[1]]&/@FixedModesG;
	For[i=1,i<=Length[FixedModes],i++,
		If[Length[FixedModes[[i]]]==4,
			AppendTo[FixedQNM,FixedModes[[i]]];AppendTo[FixedModesGPM,FixedModesG[[i]]],
			AppendTo[FixedNL,FixedModes[[i]]];AppendTo[FixedModesGNL,FixedModesG[[i]]]];
		
	];
	(*If[Length[#]==4,AppendTo[FixedQNM,#],AppendTo[FixedNL,#]]&/@FixedModes;*)
    FixedNL=Sort[#]&/@FixedNL;
    FixedModesGNL={Sort[#[[1]]],#[[2]]}&/@FixedModesGNL;
    If[#[[1,4]]==1,AppendTo[FixedModesGP,#],AppendTo[FixedModesGM,#]]&/@FixedModesGPM;
(*Find the indices for each fixed QNMs*)
	If[#[[4]]==1,AppendTo[idxFixedP,Flatten[Position[qnmp,#[[1;;3]]]]];AppendTo[FixedQNMp,#[[1;;3]]]
		,AppendTo[idxFixedM,Flatten[Position[qnmm,#[[1;;3]]]]];AppendTo[FixedQNMm,#[[1;;3]]]]&/@FixedQNM;
	(*The next two For loop is used to prepare for sorting the FixedModesG to desired order*)
	For[i=1,i<=Length[FixedModesGP],i++,
		FixedModesGP[[i]]=Join[idxFixedP[[i]],FixedModesGP[[i]]];
	];
	For[i=1,i<=Length[FixedModesGM],i++,
		FixedModesGM[[i]]=Join[idxFixedM[[i]],FixedModesGM[[i]]];
	];
	FixedModesGP=SortBy[FixedModesGP,#[[1]]&];
	FixedModesGM=SortBy[FixedModesGM,#[[1]]&];
	FixedModesGP={#[[2]],#[[3]]}&/@FixedModesGP;
	FixedModesGM={#[[2]],#[[3]]}&/@FixedModesGM;(*Sorting for qnmp and qnmm ends*)
	
	FixedQnmParamP={#,param[[2]]}&/@FixedQNMp;
	FixedQnmParamM={#,param[[2]]}&/@FixedQNMm;
	idxFixedP=Sort[Flatten[idxFixedP]];
	idxFixedM=Sort[Flatten[idxFixedM]];
(*Modify the QNM pairs so that it contains the information of a*)
	For[i=1,i<=Length[FixedQnmParamP],i++,
		qnmp[[idxFixedP[[i]]]]=FixedQnmParamP[[i]]];
	For[i=1,i<=Length[FixedQnmParamM],i++,
		qnmm[[idxFixedM[[i]]]]=FixedQnmParamM[[i]]];
	(*If there are fixed nonlinear modes*) 
	If[Length[FixedNL]!=0, 
	AppendTo[idxFixedNL,Flatten[Position[NLlist,#]]]&/@ FixedNL;
	(*Sorting the nonlinear modes in FixedModesG *)
	For[i=1,i<=Length[FixedModesGNL],i++,
		FixedModesGNL[[i]]=Join[idxFixedNL[[i]],FixedModesGNL[[i]]];
	];
	FixedModesGNL=SortBy[FixedModesGNL,#[[1]]&];
	FixedModesGNL={#[[2]],#[[3]]}&/@FixedModesGNL;
	(*Sorting ends*)
	idxFixedNL=Sort[Flatten[idxFixedNL]];
	If[#[[4]]==1,AppendTo[qnmpFixedNL,#[[1;;3]]],AppendTo[qnmmFixedNL,#[[1;;3]]]]&/@Flatten[FixedNL,1];
	qnmpFixedNL=DeleteDuplicates[qnmpFixedNL];
	qnmmFixedNL=DeleteDuplicates[qnmmFixedNL];
	AppendTo[idxFixedNlP,Flatten[Position[qnmpNL,#]]]&/@qnmpFixedNL;
	AppendTo[idxFixedNlM,Flatten[Position[qnmmNL,#]]]&/@qnmmFixedNL;
	FixedNlParamP={#,param[[2]]}&/@qnmpFixedNL;
	FixedNlParamM={#,param[[2]]}&/@qnmmFixedNL;
	For[i=1,i<=Length[FixedNlParamP],i++,
		qnmpNL[[idxFixedNlP[[i]]]]=FixedNlParamP[[i]]];
	For[i=1,i<=Length[FixedNlParamM],i++,
		qnmmNL[[idxFixedNlM[[i]]]]=FixedNlParamM[[i]]]
	];
	qnmpSet=DeleteDuplicates[Join[qnmp,qnmpNL]];
	qnmmSet=DeleteDuplicates[Join[qnmm,qnmmNL]];
	FixedModesGsorted=Join[FixedModesGP,FixedModesGM,FixedModesGNL];
	{qnmpSet,qnmmSet,{idxFixedP,idxFixedM,idxFixedNL,param[[1]]},FixedModesGsorted}                     		
]


FitInfoStruct[fit_List]:=Module[{t,\[Rho],amp,fitInfo,err2,count,singularValue,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,
	greedyInfo,qnmTemp,fixedGreedyIndex={},FixedModesG={},QQNModes,
	qnmpSet,qnmmSet,QNModespUnfixed,QNModesmUnfixed,massratioFixed,QNModespFixed,
	QNModesmFixed,QQNModesFixed,QQNModesUnfixed={},InfoTemp,i},
	(*This function is used to handle different structure of fit outputed by OverlapFit function. 
	The output of this function can be used in further data analysis. 
	The fit has a sturcture of {t,\[Rho],amp,fitInfo,err2,count}. 
	When there is quadratic modes fitted, fitInfo = {{{massratio,a,\[Theta],\[Phi]},SimModes,QNMp,QNMm,QQNModes}}.
	When the greedy algorithm is applied, fitInfo = {{massratio,a,\[Theta],\[Phi]},SimModes,QNMp,QNMm,QQNModes,greedyInfo}. 
	greedyInfo is in a form of {{{{l1,m1,n1,\[PlusMinus]1},Subscript[C, l1m1n1]},...},{mass ratio, spin magnitude}}. 
	*)
	(*Deal with the situation when fit is obtained with ReturnSingularValues=True *)
	Switch[Length[fit],
		6,
		{t,\[Rho],amp,fitInfo,err2,count}=fit,
		7,
		{t,\[Rho],amp,fitInfo,err2,count,singularValue}=fit
	];
	(* obtain the information of QNMp, QNMm, Quadratic QNM, and fixed modes. Deal with different output of OverlafFit. 
	Added for making this function compatible with qudratic QNMs and greedy algorithm.
	*)
	Switch[Length[fitInfo],
			4,
			(*FixedModesGreedy->False, NLmodesList->False*)
			{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed}=fitInfo;
			QNModesp=QNModespUnfixed;QNModesm=QNModesmUnfixed;
			qnmpSet=QNModesp;
			qnmmSet=QNModesm,
			5,
			{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,InfoTemp}=fitInfo;
			If[Length[Dimensions[InfoTemp]]==1,
			(*FixedModesGreedy->List, NLmodesList->False*)
			greedyInfo=InfoTemp;
			QNModesp=QNModespUnfixed;QNModesm=QNModesmUnfixed;
			For[i=1,i<=Length[greedyInfo[[1]]],i++,
				qnmTemp=greedyInfo[[1,i,1]];
				Switch[qnmTemp[[4]],
					1,
					AppendTo[QNModesp,qnmTemp[[1;;3]]],
					-1,
					AppendTo[QNModesm,qnmTemp[[1;;3]]]
				]
			];
			{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG}=SetGreedyModes[QNModesp,QNModesm,{},greedyInfo[[1]],greedyInfo[[2]]],
			(*FixedModesGreedy->False, NLmodesList->List*)
			QQNModesUnfixed=InfoTemp;	
			QNModesp=QNModespUnfixed;QNModesm=QNModesmUnfixed;QQNModes=QQNModesUnfixed;
			{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG}=SetGreedyModes[QNModesp,QNModesm,QQNModes,{},{massratio,a}],
			(*FixedModesGreedy->False, NLmodesList->List*)
			QQNModesUnfixed=InfoTemp;	
			QNModesp=QNModespUnfixed;QNModesm=QNModesmUnfixed;QQNModes=QQNModesUnfixed;
			{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG}=SetGreedyModes[QNModesp,QNModesm,QQNModes,{},{massratio,a}]
			]
			,
			6,
			(*FixedModesGreedy->List, NLmodesList->List*)
			{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,greedyInfo}=fitInfo;
			QNModesp=QNModespUnfixed;QNModesm=QNModesmUnfixed;QQNModes=QQNModesUnfixed;
			For[i=1,i<=Length[greedyInfo[[1]]],i++,
				qnmTemp=greedyInfo[[1,i,1]];
				Switch[Length[qnmTemp],
					4,
					Switch[qnmTemp[[4]],
					1,
					AppendTo[QNModesp,qnmTemp[[1;;3]]],
					-1,
					AppendTo[QNModesm,qnmTemp[[1;;3]]]
					],
					2,
					AppendTo[QQNModes,qnmTemp]
				]
			];
			{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG}=SetGreedyModes[QNModesp,QNModesm,QQNModes,greedyInfo[[1]],greedyInfo[[2]]]	
		];
	{{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{QNModesp,QNModesm,QQNModes},{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}}
]


(*FixRPS[rps_]:=Module[
	{amp,bhp,err,len},
	len=Length[rps[[4,3]]]+Length[rps[[4,4]]];
	amp=Transpose[Partition[rps[[3]],len],{1,5,2,3,4}];
	bhp=FlattenAt[rps[[4,1]],Position[rps[[4,1]],_?(Depth[#]==2 &),Infinity]];
	err=Transpose[Partition[rps[[5]],len],{1,5,2,3,4}];
	ReplacePart[rps,{3->amp,{4,1}->bhp,5->err}]
]*)


(* ::Section::Closed:: *)
(*End of DataRoutines Package*)


End[] (* `Private` *)


EndPackage[]
