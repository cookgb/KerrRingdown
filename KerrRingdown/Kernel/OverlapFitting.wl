(* ::Package:: *)

(* ::Chapter:: *)
(*Overlap Fitting Package*)


(* ::Section::Closed:: *)
(*Begin Overlap Fitting Package*)


BeginPackage["KerrRingdown`"]


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*Overlap Fitting*)


SetNoSpheroidalExpansion::usage=
"Calling this function causes all Overlap Fitting an Utility routines to ignore the "<>
"Spheroidal Harmonic Expansion Coefficients."


OverlapFit::usage=
"OverlapFit[BHproperties,SimModes,QNModesp,QNModesm] "<>
"Returns a list containing all of the linear-fit "<>
"information for each fit in the specified time range.\n"<>
"OverlapFit[\[Delta],\[Chi],\[Theta],\[Phi],SimModes,QNModesp,QNModesm] "<>
"Returns a list containing all of the linear-fit "<>
"information for each fit in the specified time range.\n"


FitAmplitudesTable::usage=
"FitAmplitudesTable[fit,\!\(\*SubscriptBox[\(t\), \(i\)]\)] "<>
"Creates  table of amplitudes and phases, and their standard "<>
"errors, for each QNM mode at the specified fit-start time \!\(\*SubscriptBox[\(t\), \(i\)]\)."


RestrictOverlap::usage=
"RestrictOverlap[fit] "<>
"Recompute the overlap from a fit previously "<>
"obtained by OverlapFit."


RestrictOverlap::TimeError="The specified time range is outside of the time range passed in by the OverlapFit result. 
T0=`3` cannot be smaller than `1` and TFinal=`4` cannot be greater than `2`. Please specify the option T0 and TFinal correctly."


RestrictOverlap::TimeStrideError="The following times `1` specified by RestrictOverlap do not exist in the time list passed in by the OverlapFit result. "


RelativeAmplitudes::Abort="Invalid fittime : `1`";
RelativeAmplitudes::usage=
"RelativeAmplitudes[fit,{\!\(\*SubscriptBox[\(l\), \(s\)]\),\!\(\*SubscriptBox[\(m\), \(s\)]\)},\!\(\*SubscriptBox[\(t\), \(i\)]\)] "<>
"Return the relative amplitudes contributing to signal mode "<>
"\!\(\*FormBox[SubscriptBox[\"C\", StyleBox[RowBox[{SubscriptBox[\"l\", \"s\"], SubscriptBox[\"m\", \"s\"]}],\nFontSlant->\"Italic\"]],TraditionalForm]\) "<>
"based on the the fitted coefficients at time \!\(\*SubscriptBox[\(t\), \(i\)]\).\n"<>
"RelativeAmplitudes[fit,{\!\(\*SubscriptBox[\(l\), \(s\)]\),\!\(\*SubscriptBox[\(m\), \(s\)]\)},All] "<>
"Return the relative amplitudes contributing to signal mode "<>
"\!\(\*FormBox[SubscriptBox[\(C\), \(\*SubscriptBox[\(l\), \(s\)] \*SubscriptBox[\(m\), \(s\)]\)],
TraditionalForm]\) based on the fitted coefficients at All fit times."


RemnantParameterSearch::usage=
"RemnantParameterSearch[MassRange,SpinRange,AngleRange,SimModes,QNModesp,QNModesm] "<>
"returns a structure containing the results of a set of calls to OverlapFit "<>
"with varying remnant parameters."


MaximizeOverlap::usage=
"MaximizeOverlap[\!\(\*SubscriptBox[\(\[Delta]\), \(g\)]\),\!\(\*SubscriptBox[\(\[Chi]\), \(g\)]\),\!\(\*SubscriptBox[\(\[Theta]\), \(g\)]\),\!\(\*SubscriptBox[\(t\), \(i\)]\),SimModes,QNModesp,QNModesm] "<>
"Find the optimum remnant parameters "<>
"{\!\(\*SubscriptBox[\(\[Delta]\), \(f\)]\), \!\(\*SubscriptBox[\(\[Chi]\), \(f\)]\), \!\(\*SubscriptBox[\(\[Theta]\), \(f\)]\)} to maximize "<>
"the overlap \[Rho] for fit-start time index \!\(\*SubscriptBox[\(t\), \(i\)]\).  The signal being "<>
"fit is specified by \!\(\*StyleBox[\"SimModes\", \"TI\"]\) and the QNM fitting modes are "<>
"specified by QNModesp and QNModesm.  3-dimensional maximization "<>
"is performed.\n"<>
"MaximizeOverlap[\!\(\*SubscriptBox[\(\[Delta]\), \(g\)]\),\!\(\*SubscriptBox[\(\[Chi]\), \(g\)]\),\!\(\*SubscriptBox[\(t\), \(i\)]\),SimModes,QNModesp,QNModesm] "<>
"Find the optimum remnant parameters "<>
"{\!\(\*SubscriptBox[\(\[Delta]\), \(f\)]\), \!\(\*SubscriptBox[\(\[Chi]\), \(f\)]\), 0} to maximize "<>
"the overlap \[Rho] for fit-start time index \!\(\*SubscriptBox[\(t\), \(i\)]\).  The signal being "<>
"fit is specified by \!\(\*StyleBox[\"SimModes\", \"TI\"]\) and the QNM fitting modes are "<>
"specified by QNModesp and QNModesm.  2-dimensional "<>
"maximization is performed with \[Theta] held fixed at zero"


RemnantParameterSpaceMaxOverlap::usage=
"RemnantParameterSpaceMaxOverlap[rps] "<>
"finds the maximum overlap, and associated remnant parameters, at a "<>
"sequence of fit-start times starting from the result rps obtained from "<>
"RemnantParameterSearch."


RemnantParameterSpaceMaxOverlap::TimeError="The specified time range is outside of the time range passed in by the RemnantParameterSearch result. 
T0=`3` cannot be smaller than `1` and TFinal=`4` cannot be greater than `2`. Please specify the option T0 and TFinal correctly."


RemnantParameterSpaceMaxOverlap::TimeStrideError="The following times `1` specified by RemnantParameterSpaceMaxOverlap do not exist in the time list passed in by the RemnantParameterSearch result. "


RefineMaxOverlapSequence::usage=
"RefineMaxOverlapSequence[mos,SimModes,QNModesp,QNModesm] "<>
"refines the results \!\(\*StyleBox[\"mos\", \"TI\"]\) returned by "<>
"RemnantParameterSpaceMaxOverlap "<>
"or RefineMaxOverlapSequence."


MaxOverlapSequenceAmplitudes::usage=
"MaxOverlapSequenceAmplitudes[mos,SimModes,QNModesp,QNModesm] "<>
"computes the OverlapFit using the the black-hole remnant parameter "<>
"data at every fit-start time in mos.  The QNM expansion coefficients "<>
"and their standard errors are returns and the results are presented "<>
"as amplitudes and phases."


MaxOverlapSequenceSVDInfo::usage=
"MaxOverlapSequenceSVDInfo[mos,SimModes,QNModesp,QNModesm] "<>
"recomputes the OverlapFit for each element of mos and prints "<>
"out the Singular Values."


ComputeInnerProducts::TimeError="TEnd=`2` must be greater than TFinal=`1`"
ComputeInnerProducts::TimeStrideError="Stride Length `1` must be greater than 1"
ComputeInnerProducts::TimeListError="TimeStride List is not with time range `1` to `2`"


KRFDesignMatrix::TimeStrideError="Stride Length `1` must be greater than 1"
KRFDesignMatrix::TimeListError="TimeStride List is not with time range `1` to `2`"


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Begin["`Private`"]


UseSpheroidalExpansion=True; (* Default to using Sphereoidal Expansion Coefficients *)


Protect[KRFtime,KRFC,KRF\[Omega],KRFYS,KRFYSlen,UseSpheroidalExpansion];


(* ::Section::Closed:: *)
(*Overlap Fitting*)


SetNoSpheroidalExpansion[]:=
Module[{},
	Unprotect[UseSpheroidalExpansion];
	UseSpheroidalExpansion=False;
	Protect[UseSpheroidalExpansion];
]


Options[ComputeInnerProducts]={TEnd->-1,TFinal->-2,T0->1,RestrictToSimulationSubspace->False,FitTimeStride->False,RescaleModes->False};
ComputeInnerProducts[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,OptionsPattern[]]:=
Module[{s=-2,Avec={},Belem={},Bpos={},massratio,a,\[Theta],\[Phi],t,ts,nplus,nminus,l,m,n,lp,mp,np,lpp,int,aint,rescalelist={},
		ind0=OptionValue[T0],indend=OptionValue[TEnd],indf2=OptionValue[TFinal],indf,
		subspacelpp=OptionValue[RestrictToSimulationSubspace],TimeStride=OptionValue[FitTimeStride],ListLen,TimePos,IndTime,retvec,ResModes = OptionValue[RescaleModes],tr},
	If[indend<0,indend=Length[KRFtime]+indend+1];
	If[indf2<0,indf2=Length[KRFtime]+indf2+1];
	If[indf2<indend,Null[],Message[ComputeInnerProducts::TimeError,indf2,indend];Abort[]];
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	t=Take[KRFtime,{ind0,indend}];
	tr=Take[KRFtime,{ind0,indf2}]/massratio;
	If[$MinPrecision>0,t=SetPrecision[t,$MinPrecision];tr=SetPrecision[tr,$MinPrecision];{massratio,a,\[Theta],\[Phi]}=SetPrecision[{massratio,a,\[Theta],\[Phi]},$MinPrecision]];
	ts=t/massratio;
	indf=Length[t]-(If[indend<0,Length[KRFtime]+indend+1,indend]-If[indf2<0,Length[KRFtime]+indf2+1,indf2]);
	If[Head[TimeStride]==Integer && TimeStride<1,
		Message[ComputeInnerProducts::TimeStrideError,TimeStride];Abort[]];
	If[Head[TimeStride]==List && (Min[TimeStride]<t[[1]]||Max[TimeStride]>t[[indf]]),
			Message[ComputeInnerProducts::TimeListError,t[[1]],t[[indf]]];Abort[]];
	nplus=Length[QNModesp];
	nminus=Length[QNModesm];
	Do[{l,m,n}=Qlmn; (* + modes*)
		aint=If[UseSpheroidalExpansion,
			Total[Take[KRFC@@#,{ind0,indend}]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]KRFYS[#[[1]],l,m,n]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"],
			(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
			Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"]
		];
		int=aint Exp[I Conjugate[KRF\[Omega][l,m,n]]ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[KRF\[Omega][l,m,n]]*tr]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int]
	,{Qlmn,QNModesp}];
	Do[{l,m,n}=Qlmn; (* - modes*)
		aint=If[UseSpheroidalExpansion,
			Total[Take[KRFC@@#,{ind0,indend}]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]](-1)^(l+#[[1]])KRFYS[#[[1]],l,-m,n]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"],
			(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
			Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"]
		];
		int=aint Exp[-I KRF\[Omega][l,-m,n]ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[KRF\[Omega][l,-m,n]]*tr]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int]
	,{Qlmn,QNModesm}];
	lpp=DeleteDuplicates[#[[1]]&/@SimModes]; (* restricted list of lpp values *)
	Do[{l,m,n}=QNModesp[[i]];{lp,mp,np}=QNModesp[[j]];
		If[!subspacelpp,
			(* Compute for unrestricted list of lpp values *)
			If[mp!=m,Continue[]];
			lpp=Table[l,{l,Max[Abs[s],Abs[m],Abs[mp]],Min[KRFYSlen[l,m,n]+Max[Abs[s],Abs[m]]-1,KRFYSlen[lp,mp,np]+Max[Abs[s],Abs[mp]]-1]}];
			aint=If[UseSpheroidalExpansion,
				Total[If[lp==l&&mp==m&&np==n,(Re[KRFYS[#,l,m,n]]^2+Im[KRFYS[#,l,m,n]]^2)& /@lpp,Conjugate[KRFYS[#,l,m,n]]KRFYS[#,lp,mp,np]& /@lpp],Method->"CompensatedSummation"](*The If statement needed for high precision calc*),
				(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
				If[l==lp && m==mp,Table[1,Length[ts]],Table[0,Length[ts]]]
			], 
			(* Compute for restricted list of lpp,mpp values *)
			aint=Total[Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
						KRFYS[#[[1]],lp,mp,np]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0]& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]
		];
		int=aint If[lp==l&&mp==m&&np==n,Exp[2*Im[KRF\[Omega][lp,mp,np]] ts],Exp[I(Conjugate[KRF\[Omega][l,m,n]]-KRF\[Omega][lp,mp,np]) ts]];(*The If statement needed for high precision calc*)
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]+KRF\[Omega][lp,mp,np]]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,j}];
		If[j>i,AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{j,i}]],
	{i,nplus},{j,i,nplus}];
	Do[{l,m,n}=QNModesm[[i]];{lp,mp,np}=QNModesm[[j]];
		If[!subspacelpp,
			(* Compute for unrestricted list of lpp values *)
			If[mp!=m,Continue[]];
			lpp=Table[l,{l,Max[Abs[s],Abs[m],Abs[mp]],Min[KRFYSlen[l,-m,n]+Max[Abs[s],Abs[m]]-1,KRFYSlen[lp,-mp,np]+Max[Abs[s],Abs[mp]]-1]}];
			aint=If[UseSpheroidalExpansion,
				(-1)^(l+lp)Total[If[lp==l&&mp==m&&np==n,(Re[KRFYS[#,l,-m,n]]^2+Im[KRFYS[#,l,-m,n]]^2)& /@lpp,KRFYS[#,l,-m,n]Conjugate[KRFYS[#,lp,-mp,np]]& /@lpp],Method->"CompensatedSummation"](*The If statement needed for high precision calc*),
				(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
				If[l==lp && m==mp,Table[1,Length[ts]],Table[0,Length[ts]]]
			], 
			(* Compute for restricted list of lpp,mpp values *)
			aint=Total[(-1)^(l+lp)KRFYS[#[[1]],l,-m,n]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]
						KRFYS[#[[1]],lp,-mp,np]]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0]& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]
		];
		int=aint If[lp==l&&mp==m&&np==n,Exp[2*Im[KRF\[Omega][l,-m,n]] ts],Exp[-I(KRF\[Omega][l,-m,n]-Conjugate[KRF\[Omega][lp,-mp,np]]) ts]];(*The If statement needed for high precision calc*)
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,-m,n]+KRF\[Omega][lp,-mp,np]]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+i,nplus+j}];
		If[j>i,AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+j,nplus+i}]]
	,{i,nminus},{j,i,nminus}];
	Do[{l,m,n}=QNModesp[[i]];{lp,mp,np}=QNModesm[[j]];
		If[!subspacelpp,
			(* Compute for unrestricted list of lpp values *)
			If[mp!=m,Continue[]];
			If[!subspacelpp,lpp=Table[l,{l,Max[Abs[s],Abs[m],Abs[mp]],Min[KRFYSlen[l,m,n]+Max[Abs[s],Abs[m]]-1,KRFYSlen[lp,-mp,np]+Max[Abs[s],Abs[mp]]-1]}]];
			aint=If[UseSpheroidalExpansion,
				Total[(-1)^(lp+#)Conjugate[KRFYS[#,l,m,n]KRFYS[#,lp,-mp,np]]& /@lpp,Method->"CompensatedSummation"],
				(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
				If[l==lp && m==-mp,Table[1,Length[ts]],Table[0,Length[ts]]]
			], 
			(* Compute for restricted list of lpp,mpp values *)
			aint=Total[(-1)^(#[[1]]+lp)Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]
						KRFYS[#[[1]],lp,-mp,np]]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0]& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]
		];
		int=aint Exp[I Conjugate[KRF\[Omega][l,m,n]+KRF\[Omega][lp,-mp,np]]ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]+KRF\[Omega][lp,-mp,np]]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,nplus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+j,i}],
	{i,nplus},{j,nminus}];
	int=Total[Abs[Take[KRFC @@ #,{ind0,indend}]]^2 & /@ SimModes,Method->"CompensatedSummation"];
	retvec=	{Transpose[Avec], (* Time list of A vectors *)
			Normal[SparseArray[Bpos->Flatten[#]]]&/@Transpose[Belem],(* Time list of sparse B matrices *)
			Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]], (*<Psi|Psi>*)
			Take[t,indf],
			Table[i+1,{i,indend-ind0,indend-indf2,-1}],(* Count of data/C[l,m] *)
			rescalelist};
	Switch[
	Head[TimeStride],
	Symbol,
		retvec,
	Integer,
		ListLen=Length[Transpose[Avec]];
		{retvec[[1,#]]&/@Range[1,ListLen,TimeStride], (* Time list of A vectors *)
		retvec[[2,#]]&/@Range[1,ListLen,TimeStride],(* Time list of sparse B matrices *)
		retvec[[3,#]]&/@Range[1,ListLen,TimeStride], (*<Psi|Psi>*)
		retvec[[4,#]]&/@Range[1,ListLen,TimeStride],
		retvec[[5,#]]&/@Range[1,ListLen,TimeStride],(* Count of data/C[l,m] *)
		If[ResModes,Transpose[Transpose[rescalelist][[#]]&/@Range[1,ListLen,TimeStride]],{}]},
	List,	
			(*Version 1: FitTimeStride will find the next nearest time*)
			(*IndTime=DeleteDuplicates[Flatten[SequencePosition[Take[t,indf],{x_/;x>=#},1]&/@TimeStride]];*)
			(*Version 2: FitTimeStride will find the nearest time*)
			TimePos=Flatten[(Nearest[Take[t,indf],#,1])&/@TimeStride];
			IndTime=Flatten[Position[Take[t,indf],#]&/@TimePos,2];
			{retvec[[1,#]]&/@IndTime, (* Time list of A vectors *)
			retvec[[2,#]]&/@IndTime,(* Time list of sparse B matrices *)
			retvec[[3,#]]&/@IndTime, (*<Psi|Psi>*)
			TimePos,
			retvec[[5,#]]&/@IndTime,(* Count of data/C[l,m] *)
			If[ResModes,Transpose[Transpose[rescalelist][[#]]&/@IndTime],{}]}
	]			
]


Options[KRFDesignMatrix]={TEnd->-1,TFinal->-2,T0->1,FitTimeStride->False};
KRFDesignMatrix[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,OptionsPattern[]]:=
Module[{dm,massratio,a,\[Theta],\[Phi],t,ts,i,nplus,nminus,l,m,n,
        ind0=OptionValue[T0],indend=OptionValue[TEnd],indf2=OptionValue[TFinal],indf,TimeStride=OptionValue[FitTimeStride],TimePos={},IndTime={}},
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	t=Take[KRFtime,{ind0,indend}];
	If[$MinPrecision>0,t=SetPrecision[t,$MinPrecision];{massratio,a,\[Theta],\[Phi]}=SetPrecision[{massratio,a,\[Theta],\[Phi]},$MinPrecision]];
	ts=t/massratio;
	indf=Length[t]-(If[indend<0,Length[KRFtime]+indend+1,indend]-If[indf2<0,Length[KRFtime]+indf2+1,indf2]);
	If[Head[TimeStride]==Integer && TimeStride<1,
		Message[KRFDesignMatrix::TimeStrideError,TimeStride];Abort[]];
	If[Head[TimeStride]==List && (Min[TimeStride]<t[[1]]||Max[TimeStride]>t[[indf]]),
			Message[KRFDesignMatrix::TimeListError,t[[1]],t[[indf]]];Abort[]];
	nplus=Length[QNModesp];
	nminus=Length[QNModesm];
	For[i=1,i<=nplus,++i,
		{l,m,n}=QNModesp[[i]];
		dm[i]=If[UseSpheroidalExpansion,
			Flatten[Transpose[(Exp[-I KRF\[Omega][l,m,n]ts](WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]KRFYS[#[[1]],l,m,n]))&/@SimModes]],
			(* Use this version of dm to ignore the Spheroidal Harmonic expansion coefficients *)
			Flatten[Transpose[(Exp[-I KRF\[Omega][l,m,n]ts]KroneckerDelta[l,#[[1]]])&/@SimModes]]
		];
	];
	For[i=1,i<=nminus,++i,
		{l,m,n}=QNModesm[[i]];
		dm[i+nplus]=If[UseSpheroidalExpansion,
			Flatten[Transpose[(Exp[I Conjugate[KRF\[Omega][l,-m,n]]ts](WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0](-1)^(l+#[[1]])Conjugate[KRFYS[#[[1]],l,-m,n]]))&/@SimModes]],
			(* Use this version of dm to ignore the Spheroidal Harmonic expansion coefficients *)
			Flatten[Transpose[(Exp[I Conjugate[KRF\[Omega][l,-m,n]]ts]KroneckerDelta[l,#[[1]]])&/@SimModes]]
		];
	];
	If[Head[TimeStride]==List,
	TimePos=Flatten[(Nearest[Take[t,indf],#,1])&/@TimeStride];
	IndTime=Flatten[Position[Take[t,indf],#]&/@TimePos,2]];
	{indf,Length[SimModes],Table[dm[i],{i,1,nplus+nminus}],Flatten[Transpose[Take[KRFC @@ #,{ind0,indend}]&/@SimModes]],IndTime}
]


Options[OverlapFit]=Union[{SVDWorkingPrecision->MachinePrecision,Tolerance->0,UseLeastSquares->False,ReturnSingularValues->False},Options[ComputeInnerProducts],Options[SetModeData]];
OverlapFit[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{tol=OptionValue[Tolerance],
	massratio,a,\[Theta],\[Phi]},
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	OverlapFit[massratio,a,\[Theta],\[Phi],DeleteDuplicates[SimModes],
		DeleteDuplicates[QNModesp],DeleteDuplicates[QNModesm],Evaluate@FilterRules[{opts},Options@OverlapFit]]
]


Options[OverlapFit]=Union[{SVDWorkingPrecision->MachinePrecision,Tolerance->0,UseLeastSquares->False,ReturnSingularValues->False},Options[ComputeInnerProducts],Options[SetModeData]];
OverlapFit[massratio_?NumberQ,a_?NumberQ,\[Theta]_?NumberQ,\[Phi]_?NumberQ,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{prec=OptionValue[SVDWorkingPrecision],tol=OptionValue[Tolerance],returnsv=OptionValue[ReturnSingularValues],
	Avec,Bmat,PsiPsi,t,Dmat,Qvec,count,svd,w,v,u,winv,amp,qnmp\[Delta],qnmm\[Delta],qnmp,qnmm,QNMamp,err2s,\[Rho],
	indf,nsims,dm,b,i,acol,w2,appends,nk,j,conji,dots,Rvec={},NMelem={},NMpos={},ind0,indend,ts,
	ResModes=OptionValue[RescaleModes],lp,mp,np,lm,mm,nm,rescalelist,TimeStride=OptionValue[FitTimeStride],IndTime,IndTimeDif,IndCount},
	(* Set QNM mode lists.  First Sort them. *)
	qnmp=SortBy[QNModesp,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	qnmm=SortBy[QNModesm,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	SetModeData[a,qnmp,qnmm,Evaluate@FilterRules[{opts},Options@SetModeData]];
	{Avec,Bmat,PsiPsi,t,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride]];
	Avec=SetPrecision[Avec,prec];Bmat=SetPrecision[Bmat,prec];PsiPsi=SetPrecision[PsiPsi,prec];
	Switch[OptionValue[UseLeastSquares],
	(* Eigenvalue method *)
		False,
		ind0=OptionValue[T0];
		indend=OptionValue[TFinal];
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ Bmat;
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,Avec}];	
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		\[Rho]=Sqrt[(Abs[Conjugate[#[[1]]] . (#[[3]])]^2/Abs[Conjugate[#[[3]]] . #[[2]] . (#[[3]])]&/@Transpose[{Avec,Bmat,amp}])/PsiPsi];
		If[ResModes,amp*=Transpose[rescalelist];err2s*=Transpose[rescalelist^2]],
	(* Least Squares - DesignMatix method *)
		DesignMatrix,
		{indf,nsims,dm,b,IndTime}=KRFDesignMatrix[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],FitTimeStride->OptionValue[FitTimeStride]];
		dm=Transpose[dm];
		Switch[
		Head[TimeStride],
		Symbol,
		appends=Last@Reap[For[i=1,i<=indf,++i,
			{u,w,v}=SingularValueDecomposition[SetPrecision[dm,prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[i]]]^2/Abs[Conjugate[acol] . (Bmat[[i]] . acol)]))/PsiPsi[[i]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			dm=Drop[dm,nsims]; b=Drop[b,nsims] (* remove data from latest time *)
		]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[returnsv,w:=appends[[4]]],
		Integer(* choose to fit with specific integer time stride *),
		IndCount=1;
		appends=Last@Reap[For[i=1,i<=indf,i+=TimeStride,
			{u,w,v}=SingularValueDecomposition[SetPrecision[dm,prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[IndCount]]]^2/Abs[Conjugate[acol] . (Bmat[[IndCount]] . acol)]))/PsiPsi[[IndCount]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			dm=Drop[dm,nsims*TimeStride]; b=Drop[b,nsims*TimeStride];(* remove data from latest time *)
			++IndCount 
		]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[returnsv,w:=appends[[4]]],
		List(* choose to fit with a time list *),
		dm=Drop[dm,nsims*(IndTime[[1]]-1)]; b=Drop[b,nsims*(IndTime[[1]]-1)];
		IndTimeDif=Differences[IndTime];
		IndCount=1;
		appends=Last@Reap[Do[
			{u,w,v}=SingularValueDecomposition[SetPrecision[dm,prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[IndCount]]]^2/Abs[Conjugate[acol] . (Bmat[[IndCount]] . acol)]))/PsiPsi[[IndCount]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			If[IndCount<=Length[IndTimeDif],
			dm=Drop[dm,nsims*IndTimeDif[[IndCount]]]; b=Drop[b,nsims*IndTimeDif[[IndCount]]]]; (* remove data from latest time *)
			++IndCount
		,{i,IndTime}]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[returnsv,w:=appends[[4]]]		
		],
	(* Least Squares - NormalEquation method *)
		NormalEquation,
		{indf,nsims,dm,b,IndTime}=KRFDesignMatrix[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],FitTimeStride->OptionValue[FitTimeStride]];
		nk=Dimensions[dm][[1]];
		For[i=1,i<=nk,++i,
			conji=Conjugate[dm[[i]]];
			dots=Reverse[Take[Accumulate[Reverse[conji b],Method->"CompensatedSummation"],{-(indf-1)*nsims-1,-1,nsims}]];
			Switch[Head[TimeStride],
			Symbol,
			Null,
			Integer,
			dots=dots[[#]]&/@Range[1,indf,TimeStride],
			List,
			dots=dots[[#]]&/@IndTime
			];
			AppendTo[Rvec,dots];
			For[j=i,j<=nk,++j,
				dots=Reverse[Take[Accumulate[Reverse[conji dm[[j]]],Method->"CompensatedSummation"],{-(indf-1)*nsims-1,-1,nsims}]];
				Switch[Head[TimeStride],
				Symbol,
				Null,
				Integer,
				dots=dots[[#]]&/@Range[1,indf,TimeStride],
				List,
				dots=dots[[#]]&/@IndTime
				];
				If[i==j,
				AppendTo[NMelem,Re[dots]];AppendTo[NMpos,{i,j}],
				AppendTo[NMelem,dots];AppendTo[NMpos,{i,j}];AppendTo[NMelem,Conjugate[dots]];AppendTo[NMpos,{j,i}]]
			];
		];
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ SetPrecision[Normal[SparseArray[NMpos->Flatten[#]]]&/@Transpose[NMelem],prec];
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,Transpose[Rvec]}];
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		\[Rho]=Sqrt[(Abs[Conjugate[#[[3]]] . #[[1]]]^2/Abs[Conjugate[#[[3]]] . #[[2]] . #[[3]]]&/@Transpose[{Avec,Bmat,amp}])/PsiPsi];
	];
	If[Not[returnsv],
		{t,\[Rho],SetPrecision[amp,MachinePrecision],{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm},SetPrecision[err2s,MachinePrecision],count},
		{t,\[Rho],SetPrecision[amp,MachinePrecision],{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm},SetPrecision[err2s,MachinePrecision],count,Diagonal[#]&/@w}
	]
]


Options[RestrictOverlap]=Union[{OmitModes->{{},{}}},Options[OverlapFit]];
RestrictOverlap[fit_List,opts:OptionsPattern[]]:=
Module[{prec=OptionValue[SVDWorkingPrecision],
	t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,omitp,omitm,omitposp,omitposm,omitpos,
	Avec,Bmat,PsiPsi,count,qnmp,qnmm,err2s,rescalelist,ResModes=OptionValue[RescaleModes],fittimeT0=OptionValue[T0],
	fittimeTFinal=OptionValue[TFinal],ind0,indf,tIndList,tRestrict,ampInd,ampRestrict},
	{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModesp,QNModesm},err2s,count}=fit;
	tIndList=TimeIndex[#]&/@t;
	If[fittimeT0<tIndList[[1]]||fittimeTFinal>tIndList[[-1]],Message[RestrictOverlap::TimeError,t[[1]],t[[-1]],KRFtime[[fittimeT0]],KRFtime[[fittimeTFinal]]];Abort[]];
	{omitp,omitm}=OptionValue[OmitModes];
	(* Set QNM mode lists.  First Sort them. *)
	qnmp=SortBy[QNModesp,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	qnmm=SortBy[QNModesm,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	omitposp=Table[If[MemberQ[Flatten[Position[qnmp,#]&/@omitp],i],0,1],{i,Length[qnmp]}];
	omitposm=Table[If[MemberQ[Flatten[Position[qnmm,#]&/@omitm],i],0,1],{i,Length[qnmm]}];
	omitpos=Join[omitposp,omitposm];
	SetModeData[a,qnmp,qnmm,Evaluate@FilterRules[{opts},Options@SetModeData]];
	{Avec,Bmat,PsiPsi,tRestrict,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->fittimeT0,TFinal->fittimeTFinal,TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace],RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride]];
	Avec=SetPrecision[Avec,prec];Bmat=SetPrecision[Bmat,prec];PsiPsi=SetPrecision[PsiPsi,prec];
	If[Complement[tRestrict,t]!={},Message[RestrictOverlap::TimeStrideError,Complement[tRestrict,t]];Abort[]];
	ampInd=Flatten[Position[t,#]&/@tRestrict];
	ampRestrict=amp[[#]]&/@ampInd;
	If[ResModes,ampRestrict/=Transpose[rescalelist]];
	\[Rho]=Sqrt[(Abs[Conjugate[#[[1]]] . (#[[3]]omitpos)]^2/Abs[Conjugate[#[[3]]omitpos] . #[[2]] . (#[[3]]omitpos)]&/@Transpose[{Avec,Bmat,ampRestrict}])/PsiPsi];
	{tRestrict,SetPrecision[\[Rho],MachinePrecision]}
]


Options[FitAmplitudesTable]=Options[SetModeData];
FitAmplitudesTable[fit_List,fittime_?NumberQ,opts:OptionsPattern[]]:=
Module[{t,amp,ind,ind2,QNModesp,QNModesm,qnmfixedp,qnmfixedm,massratio,a,\[Theta],\[Phi],SimModes,\[Rho],QNMamp,count,
        l,m,n,nplus,lp,mp,np,qnmp,qnmm,mode,err2s,cerr=0},
	{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModesp,QNModesm},err2s,count}=fit;
	ind=If[fittime>=t[[-1]],Length[t],If[fittime<=t[[1]],1,SequencePosition[t,{x_/;x>=fittime},1][[1,1]]]];
	ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
	SetModeData[fit[[4,1,2]],fit[[4,3]],fit[[4,4]],Evaluate@FilterRules[{opts},Options@SetModeData]];
	Do[{l,m}=Clm;
		mode=Table[0,Length[KRFtime]];
		qnmp=SpheroidalHarmonicModes[{l,m},QNModesp];
		Do[{lp,mp,np}=Qlmnp;
			mode+=If[UseSpheroidalExpansion,
				WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,Position[QNModesp,Qlmnp][[1,1]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
				(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
				KroneckerDelta[l,lp]amp[[ind,Position[QNModesp,Qlmnp][[1,1]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
			],
		{Qlmnp,qnmp}];
		nplus=Length[QNModesp];
		qnmm=SpheroidalHarmonicModes[{l,m},QNModesm];
		Do[{lp,mp,np}=Qlmnp;
			mode+=If[UseSpheroidalExpansion,
				(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,nplus+Position[QNModesm,Qlmnp][[1,1]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
				(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
				KroneckerDelta[l,lp]amp[[ind,nplus+Position[QNModesm,Qlmnp][[1,1]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
			],
		{Qlmnp,qnmm}];
		mode-=KRFC[l,m];
		mode=Take[Drop[mode,ind2-1],count[[ind]]];
		cerr+=Conjugate[mode] . mode,
	{Clm,SimModes}];
	err2s=Sqrt[err2s[[ind]]*Abs[cerr]/(2(Length[SimModes]count[[ind]]-Length[QNModesp]-Length[QNModesm]))];
	TextGrid[Join[
		{{"Mode","Amplitude","Phase/\[Pi]","\[Sigma](Amp)","\[Sigma](Phase)/\[Pi]"}},
		Flatten[#,{1,2}]&/@Transpose[{Join[{ToString[#]<>"+"}&/@QNModesp,{ToString[#]<>"-"}&/@QNModesm],
										AbsArg[Chop[#]]/{1,\[Pi]}&/@amp[[ind]],
										Transpose[{err2s,If[Im[#]!=0,1,#/\[Pi]]&/@ArcSin[err2s/Abs[amp[[ind]]]]}]}]
			],Spacings->2]
]




Options[RelativeAmplitudes]=Options[SetModeData];
RelativeAmplitudes[fit_List,lm_List,fittime_,opts:OptionsPattern[]]:=
Module[{t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,err2s,count,ind,l,m,qnmp,qnmm,lp,mp,np,j,nplus,relamps={}},
	{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModesp,QNModesm},err2s,count}=fit;
	SetModeData[a,QNModesp,QNModesm,Evaluate@FilterRules[{opts},Options@SetModeData]];
	ind=If[NumberQ[fittime],TimeIndex[fittime],
		Switch[Head[fittime],
			Symbol,If[fittime==All,Range[1,Length[t]],Message[RelativeAmplitudes::Abort,fittime];Abort[],Message[RelativeAmplitudes::Abort,fittime];Abort[]],
			_,Message[RelativeAmplitudes::Abort,fittime];Abort[]
			]
		];
	{l,m}=lm;
	qnmp=SpheroidalHarmonicModes[lm,QNModesp];
	qnmm=SpheroidalHarmonicModes[lm,QNModesm];
	If[UseSpheroidalExpansion,
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,Position[QNModesp,Qlmnp][[1,1]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] t/massratio]]],
		{Qlmnp,qnmp}];
		nplus=Length[QNModesp];
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,nplus+Position[QNModesm,Qlmnp][[1,1]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] t/massratio]]],
		{Qlmnp,qnmm}],
	(* Use this version to ignore Spheroidal Expansion expansion coefficients *)
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[KroneckerDelta[l,lp]amp[[ind,Position[QNModesp,Qlmnp][[1,1]]]]Exp[-I KRF\[Omega][lp,mp,np] t/massratio]]],
		{Qlmnp,qnmp}];
		nplus=Length[QNModesp];
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[KroneckerDelta[l,lp]amp[[ind,nplus+Position[QNModesm,Qlmnp][[1,1]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] t/massratio]]],
		{Qlmnp,qnmm}]
	];
	relamps
]


Options[RemnantParameterSearch]=Union[{},Options[OverlapFit]];
RemnantParameterSearch[MassRange_List|MassRange_?NumberQ,SpinRange_List|SpinRange_?NumberQ,AngleRange_List|AngleRange_?NumberQ,
						SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{\[Delta],\[Chi],\[Theta],sm,qnmp,qnmm,fits,tmp,mm,amp,bhp,err,cnt},
	\[Delta]=If[Head[MassRange]==List,Table[v,{v,MassRange[[1]],MassRange[[2]],MassRange[[3]]}],{MassRange},{MassRange}];
	\[Chi]=If[Head[SpinRange]==List,Table[v,{v,SpinRange[[1]],SpinRange[[2]],SpinRange[[3]]}],{SpinRange},{SpinRange}];
	\[Theta]=If[Head[AngleRange]==List,Table[v,{v,AngleRange[[1]],AngleRange[[2]],AngleRange[[3]]}],{AngleRange},{AngleRange}];
	sm=DeleteDuplicates[SimModes];
	qnmp=DeleteDuplicates[QNModesp];
	qnmm=DeleteDuplicates[QNModesm];
	fits=Array[OverlapFit[\[Delta][[#1]],\[Chi][[#2]],\[Theta][[#3]],0,sm,qnmp,qnmm,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]]&,{Length[\[Delta]],Length[\[Chi]],Length[\[Theta]]}];
	tmp=Transpose[Take[fits,All,All,All,{2}],{2,3,4,5,1}];
	mm = SetPrecision[1 - FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]],MachinePrecision];
	(*tmp=Take[fits,All,All,All,{3}];
	tmp=FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	amp = Transpose[FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]],{2,3,4,1}]; *)
	tmp=Transpose[Take[fits,All,All,All,{3}],{2,3,4,5,1,6}];
	amp = FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	tmp=Take[fits,All,All,All,{4},1];
	(*bhp = FlattenAt[tmp,Position[tmp,_?(Depth[#]==4 &),Infinity]];*)
	tmp = FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	bhp = FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	(*tmp=Take[fits,All,All,All,{5}];
	tmp=FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	err = Transpose[FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]],{2,3,4,1}];*)
	tmp=Transpose[Take[fits,All,All,All,{5}],{2,3,4,5,1,6}];
	err = FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	tmp=Transpose[Take[fits,All,All,All,{6}],{2,3,4,5,1}];
	cnt = FlattenAt[tmp,Position[tmp,_?(Depth[#]==2 &),Infinity]];
	{fits[[1,1,1,1]],mm,amp,Join[{bhp},Flatten[Take[fits,1,1,1,{4},{2,4}],4]],err,cnt}
]


(*Options[MaximizeOverlap]=Union[Complement[Options[OverlapFit],{ReturnSingularValues->False,FitTimeStride->False,T0->1,TFinal->-2}],Options[FindMaximum]];*)
MaximizeOverlap[\[Delta]g_?NumberQ,\[Chi]g_?NumberQ,\[Theta]g_?NumberQ,ti_Integer,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{sm,qnmp,qnmm,ret1=Null,ret2=Null},
	sm = DeleteDuplicates[SimModes];
	qnmp = DeleteDuplicates[QNModesp];
	qnmm = DeleteDuplicates[QNModesm];
	Off[Part::partd];
	Quiet[Check[
		ret1=SetPrecision[
			FindMaximum[OverlapFit[\[Delta]f,\[Chi]f,\[Theta]f,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],
						{{\[Delta]f,\[Delta]g,\[Delta]g+0.00001},{\[Chi]f,\[Chi]g,\[Chi]g+0.00001},{\[Theta]f,\[Theta]g,\[Theta]g+0.00001}},
						Evaluate@FilterRules[{opts},Options@FindMaximum]],
			MachinePrecision],
		Check[
		If[Head[ret1]==List,ret1[[1]]=1-ret1[[1]];
				Print["Attempt 1 failed with ",ret1],
				Print["Attempt 1 failed"],
				Print["Attempt 1 failed"]
		];
		ret2=SetPrecision[
			FindMaximum[{OverlapFit[\[Delta]f,\[Chi]f,\[Theta]f,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],0<=\[Chi]f<1},
						{{\[Delta]f,\[Delta]g},{\[Chi]f,\[Chi]g},{\[Theta]f,\[Theta]g}},Method->"Automatic",
						Evaluate@FilterRules[{opts},Options@FindMaximum]],
			MachinePrecision],
		ret2={0,{\[Delta]f->0,\[Chi]f->0,\[Theta]f->0}},
		{FindMaximum::lsbrak,FindMaximum::nrnum,Transpose::nmtx,
			SingularValueDecomposition::ovfl,General::ovfl,
			FindMaximum::cvmit}],
		{FindMaximum::lsbrak,FindMaximum::nrnum,Transpose::nmtx,
			SingularValueDecomposition::ovfl,General::ovfl,
			FindMaximum::cvmit}]];
	On[Part::partd];
	If[ret2==Null,
		Return[ret1],
		Return[ret2],
		Return[ret2]
	];
]


MaximizeOverlap[\[Delta]g_?NumberQ,\[Chi]g_?NumberQ,ti_Integer,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{sm,qnmp,qnmm,ret1=Null,ret2=Null},
	sm = DeleteDuplicates[SimModes];
	qnmp = DeleteDuplicates[QNModesp];
	qnmm = DeleteDuplicates[QNModesm];
	Off[Part::partd];
	Quiet[Check[
		ret1=SetPrecision[
			FindMaximum[OverlapFit[\[Delta]f,\[Chi]f,0,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],
						{{\[Delta]f,\[Delta]g,1.00001\[Delta]g},{\[Chi]f,\[Chi]g,1.00001\[Chi]g}},
						Evaluate@FilterRules[{opts},Options@FindMaximum]],
			MachinePrecision];
			,
		Check[
		If[Head[ret1]==List,ret1[[1]]=1-ret1[[1]];
				Print["Attempt 1 failed with ",ret1],
				Print["Attempt 1 failed"],
				Print["Attempt 1 failed"]
		];
		ret2=SetPrecision[
			FindMaximum[{OverlapFit[\[Delta]f,\[Chi]f,0,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],(0<=\[Chi]f<1)},
						{{\[Delta]f,\[Delta]g},{\[Chi]f,\[Chi]g}},Method->"Automatic",
						Evaluate@FilterRules[{opts},Options@FindMaximum]],
			MachinePrecision],
		ret2={0,{\[Delta]f->0,\[Chi]f->0}},
		{FindMaximum::lsbrak,FindMaximum::nrnum,Transpose::nmtx,
			SingularValueDecomposition::ovfl,General::ovfl,
			FindMaximum::cvmit}],
		{FindMaximum::lsbrak,FindMaximum::nrnum,Transpose::nmtx,
			SingularValueDecomposition::ovfl,General::ovfl,
			FindMaximum::cvmit}]];
	On[Part::partd];
	If[ret2==Null,
		Return[{ret1[[1]],Append[ret1[[2]],\[Theta]f->0]}],
		Return[{ret2[[1]],Append[ret2[[2]],\[Theta]f->0]}],
		Return[{ret2[[1]],Append[ret2[[2]],\[Theta]f->0]}]
	];
]


Options[RemnantParameterSpaceMaxOverlap]=Union[{FitAngle->False},Complement[Options[OverlapFit],{ReturnSingularValues->False}](*,Options[FindMaximum]*)];
RemnantParameterSpaceMaxOverlap[rps_List,opts:OptionsPattern[]]:=
Module[{Nt,params,i,pos,minparam,fit\[Theta]=OptionValue[FitAngle],maxo,ret={},tRps=rps[[1]],tIndList,
	fittimeT0=OptionValue[T0],fittimeTFinal=OptionValue[TFinal],t,tLength,TimeStride=OptionValue[FitTimeStride],
	TimePos,IndTime,IndTimeRps,tInd},
	tIndList=TimeIndex[#]&/@tRps;
	If[fittimeT0<tIndList[[1]]||fittimeTFinal>tIndList[[-1]],
		Message[RemnantParameterSpaceMaxOverlap::TimeError,tRps[[1]],tRps[[-1]],KRFtime[[fittimeT0]],KRFtime[[fittimeTFinal]]];Abort[]];
	t=Take[KRFtime,{fittimeT0,fittimeTFinal}];
	tLength=Length[t];
	Switch[
	Head[TimeStride],
	Symbol,
		t,
	Integer,
		t=t[[#]]&/@Range[1,tLength,TimeStride],
	List,	
		TimePos=Flatten[(Nearest[t,#,1])&/@TimeStride];
		IndTime=Flatten[Position[t,#]&/@TimePos,2];
		t=t[[#]]&/@IndTime		
	];
	If[Complement[t,tRps]!={},Message[RemnantParameterSpaceMaxOverlap::TimeStrideError,Complement[t,tRps]];Abort[]];
	IndTimeRps=Flatten[Position[tRps,#]&/@t];
	params=Flatten[rps[[4,1]],3];
	Nt=Length[t];
	For[i=1,i<=Nt,++i,
		tInd=IndTimeRps[[i]];
		pos=Flatten[Position[rps[[2,tInd]],Min[rps[[2,tInd]]]],1];
		(*minparam=Flatten[rps[[4,1,pos[[1]],pos[[2]],pos[[3]]]]];*)
		minparam=rps[[4,1,pos[[1]],pos[[2]],pos[[3]]]];
		maxo=If[fit\[Theta],
				MaximizeOverlap[minparam[[1]],minparam[[2]],minparam[[3]],
								TimeIndex[t[[i]]],rps[[4,2]],rps[[4,3]],rps[[4,4]],
								Evaluate@FilterRules[{opts},Options@OverlapFit]],			
				MaximizeOverlap[minparam[[1]],minparam[[2]],
								TimeIndex[t[[i]]],rps[[4,2]],rps[[4,3]],rps[[4,4]],
								Evaluate@FilterRules[{opts},Options@OverlapFit]]
			];
		Print["Time Index = ",i," : Mismatch = ",1-maxo[[1]]," \[Delta]f = ",maxo[[2,1,2]],
				" \[Chi]f = ",maxo[[2,2,2]]," \[Theta]f = ",maxo[[2,3,2]]];
		AppendTo[ret,{t[[i]],1-maxo[[1]],maxo[[2,1,2]],maxo[[2,2,2]],maxo[[2,3,2]]}];
	];
	Return[ret];
]


Options[RefineMaxOverlapSequence]=Union[{FitAngle->False,InitialGuess->None,Range->All,CheckUpdate->False},Options[OverlapFit],Options[FindMaximum]];
RefineMaxOverlapSequence[mos_List,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{Nt,i,fit\[Theta]=OptionValue[FitAngle],ig=OptionValue[InitialGuess],
		check=OptionValue[CheckUpdate],range=OptionValue[Range],
		guess,seq,rem,maxo,ret={}},
	RefineMaxOverlapSequence::InitialGuess="Invalid InitialGuess : `1`";
	{seq,rem}=TakeDrop[mos,range];
	Nt=Length[seq];
	For[i=1,i<=Nt,++i,
		Switch[ig
				,_Symbol,
				Switch[ig
					,None,guess=seq[[i,{3,4,5}]]
					,Previous,If[i==1,guess=seq[[1,{3,4,5}]],guess=ret[[i-1,{3,4,5}]]]
					,_,Message[RefineMaxOverlapSequence::InitialGuess,ig];Abort[];
					]
				,_List,guess=ig
			];		
		Print["Orig(",i,") : Time = ",seq[[i,1]]," : Mismatch = ",seq[[i,2]]," \[Delta]f = ",seq[[i,3]],
				" \[Chi]f = ",seq[[i,4]]," \[Theta]f = ",seq[[i,5]]];
		maxo=If[fit\[Theta],
				MaximizeOverlap[guess[[1]],guess[[2]],guess[[3]],
								TimeIndex[seq[[i,1]]],SimModes,QNModesp,QNModesm,
								Evaluate@FilterRules[{opts},Union[Options@OverlapFit,Options@FindMaximum]]],
				MaximizeOverlap[guess[[1]],guess[[2]],
								TimeIndex[seq[[i,1]]],SimModes,QNModesp,QNModesm,
								Evaluate@FilterRules[{opts},Union[Options@OverlapFit,Options@FindMaximum]]]
			];
		Print["New:       Time = ",KRFtime[[TimeIndex[seq[[i,1]]]]]," : Mismatch = ",1-maxo[[1]]," \[Delta]f = ",maxo[[2,1,2]],
				" \[Chi]f = ",maxo[[2,2,2]]," \[Theta]f = ",maxo[[2,3,2]]];
		If[check && 1-maxo[[1]]>=seq[[i,2]],
			Print["Discard new local solution! "];
			AppendTo[ret,seq[[i]]],
			AppendTo[ret,{seq[[i,1]],1-maxo[[1]],maxo[[2,1,2]],maxo[[2,2,2]],maxo[[2,3,2]]}]
		];
	];
	SortBy[Join[ret,rem],First]
]


MaxOverlapSequenceAmplitudes[mos_List,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{fittime,massratio,a,\[Theta],ind,ind2,fit,amp,err2,count,
		i,l,m,mode,qnmp,qnmm,j,lp,mp,np,nplus,cerr,t={},amps={},err2s={}},
	For[ind=1,ind<=Length[mos],++ind,
		fittime=mos[[ind,1]];
		ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
		{massratio,a,\[Theta]}=mos[[ind,{3,4,5}]];		fit=OverlapFit[massratio,a,\[Theta],0,DeleteDuplicates[SimModes],
			DeleteDuplicates[QNModesp],DeleteDuplicates[QNModesm],
			T0->ind2,TFinal->ind2,
			ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]];
		amp=fit[[3,1]];
		err2=fit[[5,1]];
		count=fit[[6,1]];
		cerr=0;
		For[i=1,i<=Length[SimModes],++i,
			{l,m}=SimModes[[i]];
			mode=Table[0,Length[KRFtime]];
			qnmp=SpheroidalHarmonicModes[{l,m},QNModesp];
			For[j=1,j<=Length[qnmp],++j,
				{lp,mp,np}=qnmp[[j]];
				mode+=If[UseSpheroidalExpansion,
					WignerD[{l,-m,-mp},0,\[Theta],0]amp[[Position[QNModesp,qnmp[[j]]][[1,1]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
        			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
        			KroneckerDelta[l,lp]amp[[Position[QNModesp,qnmp[[j]]][[1,1]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
        		];
			];
			nplus=Length[QNModesp];
			qnmm=SpheroidalHarmonicModes[{l,m},QNModesm];
			For[j=1,j<=Length[qnmm],++j,
				{lp,mp,np}=qnmm[[j]];
				mode+=If[UseSpheroidalExpansion,
					(-1)^(l+lp)WignerD[{l,-m,-mp},0,\[Theta],0]amp[[nplus+Position[QNModesm,qnmm[[j]]][[1,1]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
        			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
        			KroneckerDelta[l,lp]amp[[nplus+Position[QNModesm,qnmm[[j]]][[1,1]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
        		];
        	];
        	mode-=KRFC[l,m];
        	mode=Take[Drop[mode,ind2-1],count];
        	cerr+=Conjugate[mode] . mode;
        ];
        err2=Sqrt[err2*Abs[cerr]/(2(Length[SimModes]count-Length[QNModesp]-Length[QNModesm]))];
        AppendTo[t,fittime];
		AppendTo[amps,#/{1,\[Pi]}&/@AbsArg[amp]];
		AppendTo[err2s,Transpose[{err2,If[Im[#]!=0,1,#/\[Pi]]&/@ArcSin[err2/Abs[amp]]}]];
		
	];
	{t,Transpose[amps,2<->1],Transpose[err2s,2<->1]}
]


Options[MaxOverlapSequenceSVDInfo]=Union[{Range->All},Options[OverlapFit]];
MaxOverlapSequenceSVDInfo[mos_List,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{Nt,i,range=OptionValue[Range],
		guess,seq,rem,tind,of},
	{seq,rem}=TakeDrop[mos,range];
	Nt=Length[seq];
	For[i=1,i<=Nt,++i,
		tind=TimeIndex[seq[[i,1]]];
		guess=seq[[i,{3,4,5}]];
		Print["Time = ",seq[[i,1]]," : Mismatch = ",seq[[i,2]]," \[Delta]f = ",seq[[i,3]],
				" \[Chi]f = ",seq[[i,4]]," \[Theta]f = ",seq[[i,5]]];
		of=OverlapFit[guess[[1]],guess[[2]],guess[[3]],0,SimModes,QNModesp,QNModesm,
					ReturnSingularValues->True,T0->tind,TFinal->tind,
					Evaluate@FilterRules[{opts},Options@OverlapFit]];
		Print["Singular Values: ",of[[7]]];
	];
]


(* ::Section::Closed:: *)
(*End of DataRoutines Package*)


End[] (* `Private` *)


EndPackage[]
