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


SpheroidalHarmonicModes::Abort="Invalid simulation mode : `1`";
SpheroidalHarmonicModes::usage=
"SpheroidalHarmonicModes[{l,m},qnms] "<>
"Select the subset of QNMs in the list "<>
"\!\(\*StyleBox[\"qnms\", \"TI\"]\) that can overlap "<>
"with the signal mode specified by {l,m}."


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


SphericalHarmonicModes[qnm_List,sim_List]:=
Module[{s=-2},
   If[Length[qnm]!=3||qnm[[1]]<Abs[s]||Abs[qnm[[2]]]>qnm[[1]]||qnm[[3]]<0,
      Message[SphericalHarmonicModes::Abort,qnm];Abort[]];
   DeleteDuplicates[Cases[sim,x_/;x[[1]]>=Max[Abs[s],Abs[qnm[[2]]]]->x]]
]


SpheroidalHarmonicModes[sim_List,qnm_List]:=
Module[{s=-2},
   If[Length[sim]!=2||sim[[1]]<Abs[s]||Abs[sim[[2]]]>sim[[1]],
      Message[SpheroidalHarmonicModes::Abort,sim];Abort[]];
   DeleteDuplicates[Cases[Cases[qnm,x_/;x[[1]]>=Max[Abs[s],Abs[sim[[2]]]]->x],x_/;Abs[x[[2]]]<=sim[[1]]->x]]
]


Options[FitMode]={OmitModes->{{},{}}};
FitMode[fitinfo_List,lm_List,fittime_?NumberQ,OptionsPattern[]]:=
Module[{t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,qnmp,qnmm,l,m,lp,mp,np,ind,
        i,mode=Table[0,Length[KRFtime]],nplus,QNMamp,err2s,count,omitp,omitm},
	{omitp,omitm}=OptionValue[OmitModes];
	{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModesp,QNModesm},err2s,count}=fitinfo;
	{l,m}=lm;
	ind=If[fittime>=t[[-1]],Length[t],If[fittime<=t[[1]],1,SequencePosition[t,{x_/;x>=fittime},1][[1,1]]]];
	Do[{lp,mp,np}=lmn;
		If[MemberQ[omitp,lmn],Continue[]];
		mode+=If[UseSpheroidalExpansion,
			WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,Position[QNModesp,lmn][[1,1]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
			KroneckerDelta[l,lp]amp[[ind,Position[QNModesp,lmn][[1,1]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]],{lmn,SpheroidalHarmonicModes[lm,QNModesp]}];
	nplus=Length[QNModesp];
	Do[{lp,mp,np}=lmn;
		If[MemberQ[omitm,lmn],Continue[]];
		mode+=If[UseSpheroidalExpansion,
			(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,nplus+Position[QNModesm,lmn][[1,1]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
			KroneckerDelta[l,lp]amp[[ind,nplus+Position[QNModesm,lmn][[1,1]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]],{lmn,SpheroidalHarmonicModes[lm,QNModesm]}];
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
		range=OptionValue[Log10MisMatchRange],colorfunc,rng},
	{time,mismatch,amp,{BHparams,SimModes,QNModesp,QNModesm},err,count}=rps;
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


Options[OverlapSequenceAmplitudes]=Options[SetModeData];
OverlapSequenceAmplitudes[fit_List,opts:OptionsPattern[]]:=
Module[{fittime,t,\[Rho],massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,
		ind,ind2,amp,err2,count,i,l,m,mode,qnmp,qnmm,j,lp,mp,np,nplus,
		cerr,err2i,amps={},err2s={}},
	{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModesp,QNModesm},err2,count}=fit;
	SetModeData[a,QNModesp,QNModesm,Evaluate@FilterRules[{opts},Options@SetModeData]];
	For[ind=1,ind<=Length[t],++ind,
		fittime=t[[ind]];
		ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
		mode=Table[0,Length[SimModes],Length[KRFtime]];
		cerr=0;
		For[i=1,i<=Length[SimModes],++i,
			{l,m}=SimModes[[i]];
			mode=Table[0,Length[KRFtime]];
			qnmp=SpheroidalHarmonicModes[{l,m},QNModesp];
			For[j=1,j<=Length[qnmp],++j,
				{lp,mp,np}=qnmp[[j]];
				mode+=If[UseSpheroidalExpansion,
					WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,Position[QNModesp,qnmp[[j]]][[1,1]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]amp[[ind,Position[QNModesp,qnmp[[j]]][[1,1]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
				];
			];
			nplus=Length[QNModesp];
			qnmm=SpheroidalHarmonicModes[{l,m},QNModesm];
			For[j=1,j<=Length[qnmm],++j,
				{lp,mp,np}=qnmm[[j]];
				mode+=If[UseSpheroidalExpansion,
					(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,nplus+Position[QNModesm,qnmm[[j]]][[1,1]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]amp[[ind,nplus+Position[QNModesm,qnmm[[j]]][[1,1]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
				];
			];
			mode-=KRFC[l,m];
        	mode=Take[Drop[mode,ind2-1],count[[ind]]];
        	cerr+=Conjugate[mode] . mode;
        ];
        err2i=Sqrt[err2[[ind]]*Abs[cerr]/(2(Length[SimModes]count[[ind]]-Length[QNModesp]-Length[QNModesm]))];
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
