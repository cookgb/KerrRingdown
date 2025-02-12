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


Options[ComputeInnerProducts]={TEnd->-1,TFinal->-2,T0->1,RestrictToSimulationSubspace->False,FitTimeStride->False,RescaleModes->False,NLmodesList->False};
ComputeInnerProducts[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,FixedGreedyIndex_List,OptionsPattern[]]:=
Module[{s=-2,Avec={},Belem={},Bpos={},massratio,a,\[Theta],\[Phi],t,ts,nplus,nminus,l,m,n,lp,mp,np,lpp,int,aint,rescalelist={},
		ind0=OptionValue[T0],indend=OptionValue[TEnd],indf2=OptionValue[TFinal],indf,
		subspacelpp=OptionValue[RestrictToSimulationSubspace],TimeStride=OptionValue[FitTimeStride],ListLen,TimePos,IndTime,retvec,
		ResModes = OptionValue[RescaleModes],tr0,tr,NLlist=OptionValue[NLmodesList],QNModesNL,lNL,mNL,\[Omega]NL,lNLp,mNLp,\[Omega]NLp,nNL,
		fixedMassRatio,tsFixed,trFixed,tsl,trl,count,tsli,trli,tslj,trlj},
	If[indend<0,indend=Length[KRFtime]+indend+1];
	If[indf2<0,indf2=Length[KRFtime]+indf2+1];
	If[indf2<indend,Null[],Message[ComputeInnerProducts::TimeError,indf2,indend];Abort[]];
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	t=Take[KRFtime,{ind0,indend}];
	fixedMassRatio=FixedGreedyIndex[[-1]];
	tr0=Take[KRFtime,{ind0,indf2}];
	If[$MinPrecision>0,t=SetPrecision[t,$MinPrecision];tr0=SetPrecision[tr0,$MinPrecision];{massratio,a,\[Theta],\[Phi]}=SetPrecision[{massratio,a,\[Theta],\[Phi]},$MinPrecision];
	fixedMassRatio=SetPrecision[fixedMassRatio,$MinPrecision]];
	tr=tr0/massratio;
	ts=t/massratio;
	tsFixed=t/fixedMassRatio;
	trFixed=tr0/fixedMassRatio;
	indf=Length[t]-(If[indend<0,Length[KRFtime]+indend+1,indend]-If[indf2<0,Length[KRFtime]+indf2+1,indf2]);
	If[Head[TimeStride]==Integer && TimeStride<1,
		Message[ComputeInnerProducts::TimeStrideError,TimeStride];Abort[]];
	If[Head[TimeStride]==List && (Min[TimeStride]<t[[1]]||Max[TimeStride]>t[[indf]]),
			Message[ComputeInnerProducts::TimeListError,t[[1]],t[[indf]]];Abort[]];
	(*Given a NLlist which component has the form of {{l1,m1,n1,\[PlusMinus]},{l2,m2,n2,\[PlusMinus]}},
	the code below finds the frequency of the quadratic mode. *)	
	If[Head[NLlist]==List,
		QNModesNL=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@NLlist];
		nNL=Length[QNModesNL]
		];
	nplus=Length[QNModesp];
	nminus=Length[QNModesm];
	
	count=1;
	Do[
		If[MemberQ[FixedGreedyIndex[[1]],count],
			tsl:=tsFixed;trl:=trFixed,
			tsl:=ts;trl:=tr];
		{l,m,n}=Qlmn; (* + modes*)
		aint=If[UseSpheroidalExpansion,
			Total[Take[KRFC@@#,{ind0,indend}]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]KRFYS[#[[1]],l,m,n]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"],
			(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
			Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"]
		];
		int=aint Exp[I Conjugate[KRF\[Omega][l,m,n]]tsl];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[KRF\[Omega][l,m,n]]*trl]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int];
		count++
	,{Qlmn,QNModesp}];
	
	count=1;
	Do[
		If[MemberQ[FixedGreedyIndex[[2]],count],
			tsl:=tsFixed;trl:=trFixed,
			tsl:=ts;trl:=tr];
		{l,m,n}=Qlmn; (* - modes*)
		aint=If[UseSpheroidalExpansion,
			Total[Take[KRFC@@#,{ind0,indend}]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]](-1)^(l+#[[1]])KRFYS[#[[1]],l,-m,n]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"],
			(* Use this version of aint to ignore the Spheroidal Harmonic expansion coefficients *)
			Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]&/@SphericalHarmonicModes[Qlmn,SimModes],Method->"CompensatedSummation"]
		];
		int=aint Exp[-I KRF\[Omega][l,-m,n]tsl];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[KRF\[Omega][l,-m,n]]*trl]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int];
		count++
	,{Qlmn,QNModesm}];
	
	count=1;
	If[Head[NLlist]==List,
	Do[
		If[MemberQ[FixedGreedyIndex[[3]],count],
			tsl:=tsFixed;trl:=trFixed,
			tsl:=ts;trl:=tr];
		{l,m,\[Omega]NL}=Qlmn;(* nonlinear modes*)
		aint=Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]&/@SphericalHarmonicModesNL[Qlmn,SimModes],Method->"CompensatedSummation"];
		int=aint Exp[I Conjugate[\[Omega]NL]tsl];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[\[Omega]NL]*trl]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int];
		count++
	,{Qlmn,QNModesNL}]];
	
	
	lpp=DeleteDuplicates[#[[1]]&/@SimModes]; (* restricted list of lpp values *)
	Do[
	If[MemberQ[FixedGreedyIndex[[1]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
	If[MemberQ[FixedGreedyIndex[[1]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
	{l,m,n}=QNModesp[[i]];{lp,mp,np}=QNModesp[[j]];
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
			aint=Total[If[lp==l&&mp==m&&np==n,
							((Re[KRFYS[#[[1]],l,m,n]]^2+Im[KRFYS[#[[1]],l,m,n]]^2)Abs[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]^2)& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],
							(Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KRFYS[#[[1]],lp,mp,np]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&]],Method->"CompensatedSummation"]
		];
		int=aint If[lp==l&&mp==m&&np==n,Exp[2*Im[KRF\[Omega][lp,mp,np]] tsli],Exp[I(Conjugate[KRF\[Omega][l,m,n]]*tsli-KRF\[Omega][lp,mp,np]*tslj)]];(*The If statement needed for high precision calc*)
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]]*trli-Im[KRF\[Omega][lp,mp,np]]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,j}];
		If[j>i,AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{j,i}]],
	{i,nplus},{j,i,nplus}];
	
	Do[
		If[MemberQ[FixedGreedyIndex[[2]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
		If[MemberQ[FixedGreedyIndex[[2]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
		{l,m,n}=QNModesm[[i]];{lp,mp,np}=QNModesm[[j]];
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
			aint=Total[If[lp==l&&mp==m&&np==n,
							((-1)^(l+lp)(Re[KRFYS[#[[1]],l,-m,n]]^2+Im[KRFYS[#[[1]],l,-m,n]]^2)Abs[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]^2)& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],
							((-1)^(l+lp)KRFYS[#[[1]],l,-m,n]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]
							KRFYS[#[[1]],lp,-mp,np]]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&]],Method->"CompensatedSummation"]
							
		];
		int=aint If[lp==l&&mp==m&&np==n,Exp[2*Im[KRF\[Omega][l,-m,n]] tsli],Exp[-I(KRF\[Omega][l,-m,n]*tsli-Conjugate[KRF\[Omega][lp,-mp,np]]*tslj)]];(*The If statement needed for high precision calc*)
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,-m,n]]*trli-Im[KRF\[Omega][lp,-mp,np]]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+i,nplus+j}];
		If[j>i,AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+j,nplus+i}]]
	,{i,nminus},{j,i,nminus}];	
	
	Do[
		If[MemberQ[FixedGreedyIndex[[1]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
		If[MemberQ[FixedGreedyIndex[[2]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
		{l,m,n}=QNModesp[[i]];{lp,mp,np}=QNModesm[[j]];
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
		int=aint Exp[I (Conjugate[KRF\[Omega][l,m,n]]*tsli+Conjugate[KRF\[Omega][lp,-mp,np]]*tslj)];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]]*trli-Im[KRF\[Omega][lp,-mp,np]]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,nplus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+j,i}],
	{i,nplus},{j,nminus}];	
	If[Head[NLlist]==List,
	Do[
		If[MemberQ[FixedGreedyIndex[[1]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
		If[MemberQ[FixedGreedyIndex[[3]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
		{l,m,n}=QNModesp[[i]];{lNL,mNL,\[Omega]NL}=QNModesNL[[j]];(*+modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=m,Continue[]];
			aint=If[UseSpheroidalExpansion,Conjugate[KRFYS[lNL,l,m,n]],KroneckerDelta[lNL,l]],
			aint=Total[(Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]				
			];
		int=aint Exp[I(Conjugate[KRF\[Omega][l,m,n]]*tsli-\[Omega]NL*tslj)];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]]*trli-Im[\[Omega]NL]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,i}],
	{i,nplus},{j,nNL}];
	Do[
		If[MemberQ[FixedGreedyIndex[[2]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
		If[MemberQ[FixedGreedyIndex[[3]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
		{l,m,n}=QNModesm[[i]];{lNL,mNL,\[Omega]NL}=QNModesNL[[j]];(*-modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=m,Continue[]];
			aint=If[UseSpheroidalExpansion,(-1)^(l+lNL)KRFYS[lNL,l,-m,n],KroneckerDelta[lNL,l]],
			aint=Total[((-1)^(l+lp)KRFYS[#[[1]],l,-m,n]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]];
		int=aint Exp[-I(KRF\[Omega][l,-m,n]*tsli+\[Omega]NL*tslj)];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,-m,n]]*trli-Im[\[Omega]NL]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,nplus+i}],
	{i,nminus},{j,nNL}];
	Do[
		If[MemberQ[FixedGreedyIndex[[3]],i],tsli:=tsFixed;trli:=trFixed,tsli:=ts;trli:=tr];
		If[MemberQ[FixedGreedyIndex[[3]],j],tslj:=tsFixed;trlj:=trFixed,tslj:=ts;trlj:=tr];
		{lNL,mNL,\[Omega]NL}=QNModesNL[[i]];{lNLp,mNLp,\[Omega]NLp}=QNModesNL[[j]];(*nonlinear modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=mNLp && lNL!=lNLp,Continue[]];
			aint=1,
			aint=Total[(Conjugate[WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] KroneckerDelta[#[[1]],lNLp] WignerD[{#[[1]],-#[[2]],-mNLp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]
		];	
		int=aint Exp[I(Conjugate[\[Omega]NL]*tsli-\[Omega]NLp*tslj)];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[\[Omega]NL]*trli-Im[\[Omega]NLp]*trlj]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+nminus+i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,nplus+nminus+i}],
	{i,nNL},{j,nNL}]];	
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


Options[ComputeInnerProducts]={TEnd->-1,TFinal->-2,T0->1,RestrictToSimulationSubspace->False,FitTimeStride->False,RescaleModes->False,NLmodesList->False};
ComputeInnerProducts[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,OptionsPattern[]]:=
Module[{s=-2,Avec={},Belem={},Bpos={},massratio,a,\[Theta],\[Phi],t,ts,nplus,nminus,l,m,n,lp,mp,np,lpp,int,aint,rescalelist={},
		ind0=OptionValue[T0],indend=OptionValue[TEnd],indf2=OptionValue[TFinal],indf,
		subspacelpp=OptionValue[RestrictToSimulationSubspace],TimeStride=OptionValue[FitTimeStride],ListLen,TimePos,IndTime,retvec,
		ResModes = OptionValue[RescaleModes],tr0,tr,NLlist=OptionValue[NLmodesList],QNModesNL,lNL,mNL,\[Omega]NL,lNLp,mNLp,\[Omega]NLp,nNL},
	If[indend<0,indend=Length[KRFtime]+indend+1];
	If[indf2<0,indf2=Length[KRFtime]+indf2+1];
	If[indf2<indend,Null[],Message[ComputeInnerProducts::TimeError,indf2,indend];Abort[]];
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	t=Take[KRFtime,{ind0,indend}];
	tr0=Take[KRFtime,{ind0,indf2}];
	If[$MinPrecision>0,t=SetPrecision[t,$MinPrecision];tr0=SetPrecision[tr0,$MinPrecision];{massratio,a,\[Theta],\[Phi]}=SetPrecision[{massratio,a,\[Theta],\[Phi]},$MinPrecision]];
	tr=tr0/massratio;
	ts=t/massratio;
	indf=Length[t]-(If[indend<0,Length[KRFtime]+indend+1,indend]-If[indf2<0,Length[KRFtime]+indf2+1,indf2]);
	If[Head[TimeStride]==Integer && TimeStride<1,
		Message[ComputeInnerProducts::TimeStrideError,TimeStride];Abort[]];
	If[Head[TimeStride]==List && (Min[TimeStride]<t[[1]]||Max[TimeStride]>t[[indf]]),
			Message[ComputeInnerProducts::TimeListError,t[[1]],t[[indf]]];Abort[]];
	(*Given a NLlist which component has the form of {{l1,m1,n1,\[PlusMinus]},{l2,m2,n2,\[PlusMinus]}},
	the code below finds the frequency of the quadratic mode. *)	
	If[Head[NLlist]==List,
		QNModesNL=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@NLlist];
		nNL=Length[QNModesNL]
		];
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
	If[Head[NLlist]==List,
	Do[{l,m,\[Omega]NL}=Qlmn;(* nonlinear modes*)
		aint=Total[Take[KRFC@@#,{ind0,indend}]KroneckerDelta[l,#[[1]]]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]&/@SphericalHarmonicModesNL[Qlmn,SimModes],Method->"CompensatedSummation"];
		int=aint Exp[I Conjugate[\[Omega]NL]ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes,AppendTo[rescalelist,Exp[-Im[\[Omega]NL]*tr]]; int*=rescalelist[[-1]]];
		AppendTo[Avec,int]
	,{Qlmn,QNModesNL}]];
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
			aint=Total[If[lp==l&&mp==m&&np==n,
							((Re[KRFYS[#[[1]],l,m,n]]^2+Im[KRFYS[#[[1]],l,m,n]]^2)Abs[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]^2)& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],
							(Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KRFYS[#[[1]],lp,mp,np]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&]],Method->"CompensatedSummation"]
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
			aint=Total[If[lp==l&&mp==m&&np==n,
							((-1)^(l+lp)(Re[KRFYS[#[[1]],l,-m,n]]^2+Im[KRFYS[#[[1]],l,-m,n]]^2)Abs[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]^2)& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],
							((-1)^(l+lp)KRFYS[#[[1]],l,-m,n]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]
							KRFYS[#[[1]],lp,-mp,np]]WignerD[{#[[1]],-#[[2]],-mp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&]],Method->"CompensatedSummation"]
							
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
	If[Head[NLlist]==List,
	Do[{l,m,n}=QNModesp[[i]];{lNL,mNL,\[Omega]NL}=QNModesNL[[j]];(*+modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=m,Continue[]];
			aint=If[UseSpheroidalExpansion,Conjugate[KRFYS[lNL,l,m,n]],KroneckerDelta[lNL,l]],
			aint=Total[(Conjugate[KRFYS[#[[1]],l,m,n]WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]				
			];
		int=aint Exp[I(Conjugate[KRF\[Omega][l,m,n]]-\[Omega]NL) ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,m,n]+\[Omega]NL]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,i}],
	{i,nplus},{j,nNL}];
	Do[{l,m,n}=QNModesm[[i]];{lNL,mNL,\[Omega]NL}=QNModesNL[[j]];(*-modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=m,Continue[]];
			aint=If[UseSpheroidalExpansion,(-1)^(l+lNL)KRFYS[lNL,l,-m,n],KroneckerDelta[lNL,l]],
			aint=Total[((-1)^(l+lp)KRFYS[#[[1]],l,-m,n]Conjugate[WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]];
		int=aint Exp[-I(KRF\[Omega][l,-m,n]+\[Omega]NL) ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[KRF\[Omega][l,-m,n]+\[Omega]NL]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,nplus+i}],
	{i,nminus},{j,nNL}];
	Do[{lNL,mNL,\[Omega]NL}=QNModesNL[[i]];{lNLp,mNLp,\[Omega]NLp}=QNModesNL[[j]];(*nonlinear modes and nonlinear modes*)
		If[!subspacelpp,
			If[mNL!=mNLp && lNL!=lNLp,Continue[]];
			aint=1,
			aint=Total[(Conjugate[WignerD[{#[[1]],-#[[2]],-mNL},\[Phi],\[Theta],0]]
							KroneckerDelta[#[[1]],lNL] KroneckerDelta[#[[1]],lNLp] WignerD[{#[[1]],-#[[2]],-mNLp},\[Phi],\[Theta],0])& /@Select[SimModes,#[[1]]>=Max[Abs[m],Abs[mp]]&],Method->"CompensatedSummation"]
		];	
		int=aint Exp[I(Conjugate[\[Omega]NL]-\[Omega]NLp) ts];
		int=Reverse[Take[Accumulate[Reverse[(Drop[int,1]+Drop[int,-1])/2],Method->"CompensatedSummation"],-indf]];
		If[ResModes, int *=  Exp[-Im[\[Omega]NL+\[Omega]NLp]*tr]];
		AppendTo[Belem,int];AppendTo[Bpos,{nplus+nminus+i,nplus+nminus+j}];
		AppendTo[Belem,Conjugate[int]];AppendTo[Bpos,{nplus+nminus+j,nplus+nminus+i}],
	{i,nNL},{j,nNL}]];	
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


Options[KRFDesignMatrix]={TEnd->-1,TFinal->-2,T0->1,FitTimeStride->False,NLmodesList->False};
KRFDesignMatrix[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,OptionsPattern[]]:=
Module[{dm,massratio,a,\[Theta],\[Phi],t,ts,i,nplus,nminus,l,m,n,
        ind0=OptionValue[T0],indend=OptionValue[TEnd],indf2=OptionValue[TFinal],indf,TimeStride=OptionValue[FitTimeStride],
        TimePos={},IndTime={},NLlist=OptionValue[NLmodesList],QNModesNL,nNL,lNL,mNL,\[Omega]NL,\[Omega]Info={}},
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	t=Take[KRFtime,{ind0,indend}];
	If[$MinPrecision>0,t=SetPrecision[t,$MinPrecision];{massratio,a,\[Theta],\[Phi]}=SetPrecision[{massratio,a,\[Theta],\[Phi]},$MinPrecision]];
	ts=t/massratio;
	indf=Length[t]-(If[indend<0,Length[KRFtime]+indend+1,indend]-If[indf2<0,Length[KRFtime]+indf2+1,indf2]);
	If[Head[TimeStride]==Integer && TimeStride<1,
		Message[KRFDesignMatrix::TimeStrideError,TimeStride];Abort[]];
	If[Head[TimeStride]==List && (Min[TimeStride]<t[[1]]||Max[TimeStride]>t[[indf]]),
			Message[KRFDesignMatrix::TimeListError,t[[1]],t[[indf]]];Abort[]];
	If[Head[NLlist]==List,
		QNModesNL=DeleteDuplicates[{#[[1,1]]+#[[2,1]],#[[1,2]]+#[[2,2]],
				If[#[[1,4]]==1,
					If[#[[2,4]]==1,
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						KRF\[Omega][#[[1,1]],#[[1,2]],#[[1,3]]]-Conjugate[KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					],
					If[#[[2,4]]==1,
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]]+KRF\[Omega][#[[2,1]],#[[2,2]],#[[2,3]]],
						-Conjugate[KRF\[Omega][#[[1,1]],-#[[1,2]],#[[1,3]]]+KRF\[Omega][#[[2,1]],-#[[2,2]],#[[2,3]]]]
					]
					]}&/@NLlist];
		nNL=Length[QNModesNL]
	];
	nplus=Length[QNModesp];
	nminus=Length[QNModesm];
	For[i=1,i<=nplus,++i,
		{l,m,n}=QNModesp[[i]];
		dm[i]=If[UseSpheroidalExpansion,
			Flatten[Transpose[(Exp[-I KRF\[Omega][l,m,n]ts](WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]KRFYS[#[[1]],l,m,n]))&/@SimModes]],
			(* Use this version of dm to ignore the Spheroidal Harmonic expansion coefficients *)
			Flatten[Transpose[(Exp[-I KRF\[Omega][l,m,n]ts]KroneckerDelta[l,#[[1]]])&/@SimModes]]
		];
		AppendTo[\[Omega]Info,KRF\[Omega][l,m,n]]
	];
	For[i=1,i<=nminus,++i,
		{l,m,n}=QNModesm[[i]];
		dm[i+nplus]=If[UseSpheroidalExpansion,
			Flatten[Transpose[(Exp[I Conjugate[KRF\[Omega][l,-m,n]]ts](WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0](-1)^(l+#[[1]])Conjugate[KRFYS[#[[1]],l,-m,n]]))&/@SimModes]],
			(* Use this version of dm to ignore the Spheroidal Harmonic expansion coefficients *)
			Flatten[Transpose[(Exp[I Conjugate[KRF\[Omega][l,-m,n]]ts]KroneckerDelta[l,#[[1]]])&/@SimModes]]
		];
		AppendTo[\[Omega]Info,-Conjugate[KRF\[Omega][l,-m,n]]]
	];
	If[Head[NLlist]==List,
	For[i=1,i<=nNL,++i, (*nonlinear modes*)
		{lNL,mNL,\[Omega]NL}=QNModesNL[[i]];
		dm[i+nplus+nminus]=If[UseSpheroidalExpansion,
			Flatten[Transpose[(Exp[-I \[Omega]NL ts](WignerD[{#[[1]],-#[[2]],-m},\[Phi],\[Theta],0]KroneckerDelta[lNL,#[[1]]]))&/@SimModes]],
			(* Use this version of dm to ignore the Spheroidal Harmonic expansion coefficients *)
			Flatten[Transpose[(Exp[-I \[Omega]NL ts]KroneckerDelta[lNL,#[[1]]])&/@SimModes]]
		];
		AppendTo[\[Omega]Info,\[Omega]NL]
	]];
	If[Head[TimeStride]==List,
	TimePos=Flatten[(Nearest[Take[t,indf],#,1])&/@TimeStride];
	IndTime=Flatten[Position[Take[t,indf],#]&/@TimePos,2]];
	{indf,Length[SimModes],If[Head[NLlist]==List,Table[dm[i],{i,1,nplus+nminus+nNL}],Table[dm[i],{i,1,nplus+nminus}],Table[dm[i],{i,1,nplus+nminus}]],Flatten[Transpose[Take[KRFC @@ #,{ind0,indend}]&/@SimModes]],IndTime,\[Omega]Info}
]


Options[OverlapFit]=Union[{SVDWorkingPrecision->MachinePrecision,Tolerance->0,UseLeastSquares->False,ReturnSingularValues->False,FixedModesGreedy->False,FullMismatch->False},Options[ComputeInnerProducts],Options[SetModeData]];
OverlapFit[BHproperties_List,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{tol=OptionValue[Tolerance],
	massratio,a,\[Theta],\[Phi]},
	{massratio,a,\[Theta],\[Phi]}=BHproperties;
	OverlapFit[massratio,a,\[Theta],\[Phi],DeleteDuplicates[SimModes],
		DeleteDuplicates[QNModesp],DeleteDuplicates[QNModesm],Evaluate@FilterRules[{opts},Options@OverlapFit]]
]


Options[OverlapFit]=Union[{SVDWorkingPrecision->MachinePrecision,Tolerance->0,UseLeastSquares->False,ReturnSingularValues->False,FixedModesGreedy->False,FullMismatch->False},Options[ComputeInnerProducts],Options[SetModeData]];
(* When the greedy algorithm is used, the info of QNMs to be fixed is passed in as option FixedModesGreedy, 
which is in a form of {{{{l1,m1,n1,\[PlusMinus]},Subscript[C, l1m1n1]},...},{mass ratio, spin magnitude}}. 
The remnant parameters are used to fix their QNM frequencies. *)
OverlapFit[massratio_?NumberQ,a_?NumberQ,\[Theta]_?NumberQ,\[Phi]_?NumberQ,SimModes_List,QNModesp_List,QNModesm_List,opts:OptionsPattern[]]:=
Module[{prec=OptionValue[SVDWorkingPrecision],tol=OptionValue[Tolerance],returnsv=OptionValue[ReturnSingularValues],
	Avec,Bmat,PsiPsi,t,Dmat,Qvec,count,svd,w,v,u,winv,amp,qnmp\[Delta],qnmm\[Delta],qnmp,qnmm,QNMamp,err2s,\[Rho]
	,nsims,dm,b,i,acol,w2,appends,nk,j,conji,dots,Rvec={},NMelem={},NMpos={},
	indf,ind0=OptionValue[T0],indf2=OptionValue[TFinal],tr,ts,
	ResModes=OptionValue[RescaleModes],lp,mp,np,lm,mm,nm,rescalelist,TimeStride=OptionValue[FitTimeStride],IndTime,IndTimeDif,IndCount,
	FixedModesG,fixedModesGreedy=OptionValue[FixedModesGreedy],FixedModes={},FixedCoef,
	idxFixed={},idxUnknown,nplus,nminus,
	AvecK,AvecU,Btrans,BmatK,BKelem={},BKpos={},BmatU,BUelem={},BUpos={},BmatC,BCelem={},BCpos={},\[Alpha]mat,
	BmatCconj,nUnknown,nFixed,PsiPsiU,qnmpInd={},qnmmInd={},FixedCount,UnknownCount,ampTotal={},ampFixed,ampUnknown,
	rescalelistUnknown,rescalelistFixed,\[Rho]T,(*after that is the variables used for nonlinear modes*)
	NLlist=OptionValue[NLmodesList],qnmpNL={},qnmmNL={},qnmpSet,qnmmSet,FixedModesNL={},NLInd={},
	paramGreedy,fixedGreedyIndex,Binv,fullMismatch=OptionValue[FullMismatch],\[Rho]2,\[Omega]Info,dmCurr},
	If[indf2<0,indf2=Length[KRFtime]+indf2+1]; 
	tr=Take[KRFtime,{ind0,indf2}]/massratio;(*This tr is used for rescaling in least-square*)
	(* Set QNM mode lists.  First Sort them. *)
	(*make sure a list or false. Make sure 4 index is +1 or -1*)
	If[Head[NLlist]==List||NLlist==False,Null,
      Message[OverlapFit::NonlinearModesOptionError,NLlist];Abort[]];
    If[Head[NLlist]==List,If[#[[4]]!=1&&#[[4]]!=-1,Message[OverlapFit::NonlinearModesListError];Abort[]]&/@Flatten[NLlist,1]];
    If[Head[NLlist]==List,NLlist=DeleteDuplicates[Sort[#]&/@NLlist]];(*Sort each nonlinear modes list to make Position function works better*)
	If[Head[NLlist]==List,If[#[[4]]==1,AppendTo[qnmpNL,#[[1;;3]]],AppendTo[qnmmNL,#[[1;;3]]]]&/@DeleteDuplicates[Flatten[NLlist,1]]];
	qnmp\[Delta]=SortBy[QNModesp,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	qnmm\[Delta]=SortBy[QNModesm,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	If[Head[fixedModesGreedy]==List,
	(*If there is a fourth element in qnmpSet, just send an error message. It does not allow
		testing over fitting right now.*)
		paramGreedy=fixedModesGreedy[[-1]];
		FixedModesG=fixedModesGreedy[[1]];
		{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG}=SetGreedyModes[qnmp\[Delta],qnmm\[Delta],NLlist,FixedModesG,paramGreedy],
		(*SetGreedyModes can also give out the indicies of fixed modes (fixedGreedyIndex). *)
		qnmpSet=DeleteDuplicates[Join[qnmp\[Delta],qnmpNL]];
		qnmmSet=DeleteDuplicates[Join[qnmm\[Delta],qnmmNL]],
		qnmpSet=DeleteDuplicates[Join[qnmp\[Delta],qnmpNL]];
		qnmmSet=DeleteDuplicates[Join[qnmm\[Delta],qnmmNL]]	
	];
	SetModeData[a,qnmpSet,qnmmSet,Evaluate@FilterRules[{opts},Options@SetModeData]];	
	{qnmp,qnmm} = {(Take[#,3]&/@qnmp\[Delta]),(Take[#,3]&/@qnmm\[Delta])};
	nplus=Length[qnmp];
	nminus=Length[qnmm];
	If[Head[fixedModesGreedy]==List,
	{Avec,Bmat,PsiPsi,t,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,fixedGreedyIndex,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->OptionValue[NLmodesList]],
	{Avec,Bmat,PsiPsi,t,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->OptionValue[NLmodesList]],
	{Avec,Bmat,PsiPsi,t,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->OptionValue[NLmodesList]]	
	];
		
	Avec=SetPrecision[Avec,prec];Bmat=SetPrecision[Bmat,prec];PsiPsi=SetPrecision[PsiPsi,prec];
	rescalelist=SetPrecision[rescalelist,prec];
	
	(*If the Greedy Algorithm is used*) 
	If[Head[fixedModesGreedy]==List,
	If[Head[NLlist]==List(*If there are nonlinear modes in the list*), 
		(*obtain the indices of QNMs and nonlinear modes *)
		idxUnknown=Range[1,nplus+nminus+Length[NLlist]];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		NLInd={#}&/@fixedGreedyIndex[[3]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus,fixedGreedyIndex[[3]]+nplus+nminus];
		idxFixed={#}&/@idxFixed;
		NLlist = Delete[NLlist,NLInd],
		(*If no nonlinear modes in the list*)
		idxUnknown=Range[1,nplus+nminus];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus];
		idxFixed={#}&/@idxFixed,
		idxUnknown=Range[1,nplus+nminus];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus];
		idxFixed={#}&/@idxFixed
	];
	qnmp = Delete[qnmp,qnmpInd];
	qnmm = Delete[qnmm,qnmmInd];
	idxUnknown=Delete[idxUnknown,idxFixed];
	idxFixed=Flatten[idxFixed];
	nUnknown=Length[idxUnknown];
	nFixed=Length[idxFixed];	
	(*Seperate unknown and known parts from the A vector *)
	AvecK=Transpose[Part[Transpose[Avec],idxFixed]];
	AvecU=Transpose[Part[Transpose[Avec],idxUnknown]];
	(*Seperate B fixed from B matrix*)
	Btrans=Transpose[Bmat,{3,1,2}];
	Do[AppendTo[BKelem,Btrans[[idxFixed[[i]],idxFixed[[j]]]]];AppendTo[BKpos,{i,j}];
		If[j>i,AppendTo[BKelem,Conjugate[Btrans[[idxFixed[[i]],idxFixed[[j]]]]]];AppendTo[BKpos,{j,i}]]
		,{i,nFixed},{j,i,nFixed}];
	BmatK=Normal[SparseArray[BKpos->Flatten[#]]]&/@Transpose[BKelem];
	(*Seperate B unknown from B matrix *)
	Do[AppendTo[BUelem,Btrans[[idxUnknown[[i]],idxUnknown[[j]]]]];AppendTo[BUpos,{i,j}];
		If[j>i,AppendTo[BUelem,Conjugate[Btrans[[idxUnknown[[i]],idxUnknown[[j]]]]]];AppendTo[BUpos,{j,i}]]
		,{i,nUnknown},{j,i,nUnknown}];
	BmatU=Normal[SparseArray[BUpos->Flatten[#]]]&/@Transpose[BUelem];	
	(*Seperate B fixed and unknown overlap part from B matrix *)
	Do[AppendTo[BCelem,Btrans[[idxFixed[[i]],idxUnknown[[j]]]]];AppendTo[BCpos,{i,j}]
		,{i,nFixed},{j,nUnknown}];
	BmatC=Normal[SparseArray[BCpos->Flatten[#]]]&/@Transpose[BCelem];
	BmatCconj=ConjugateTranspose[#]&/@BmatC;
	(*Set up fixed Coefficients. *)
	FixedCoef=Transpose[(Table[1,{i,Length[Avec]}]*#[[2]])&/@FixedModesG];
	(* Get the fixed and unknown parts of the rescalelist. If rescaling is applied, the 
	fixed coefficients should also be changed to match the rescaling *)
	If[ResModes,
		rescalelistUnknown=Part[rescalelist,idxUnknown];
		rescalelistFixed=Part[rescalelist,idxFixed];
		FixedCoef/=Transpose[rescalelistFixed]];
	(*Calculate the unknown part of PsiPsi*)
	\[Alpha]mat=(#[[1]]-#[[2]] . #[[3]])&/@Transpose[{AvecU,BmatCconj,FixedCoef}];
	PsiPsiU=(#[[1]]-2*Re[Conjugate[#[[4]]] . #[[3]]]+Conjugate[#[[4]]] . (#[[2]] . #[[4]]))&/@Transpose[{PsiPsi,BmatK,AvecK,FixedCoef}]
	];
	
	(*Note that the Greedy Algorithm is applied to Eigenvalue Method.*)
	Switch[OptionValue[UseLeastSquares],
	(* Eigenvalue method *)
		False,
		If[Head[FixedModesG]!=List,
		(*If Greedy Algorithm is not used*)
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ Bmat;
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,Avec}];
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		(*\[Rho]=Sqrt[Abs[((Conjugate[#[[1]]] . #[[2]])&/@Transpose[{Avec,amp}])/PsiPsi]];
		Both forms of \[Rho] are equivalent with each other theoretically. 
		The overlap function used below can give overlaps with lower numerical noise, 
		which is very useful for nonlinear fitting scenario. 
		*)
		\[Rho]=Sqrt[(Abs[Conjugate[#[[1]]] . (#[[3]])]^2/Abs[Conjugate[#[[3]]] . #[[2]] . (#[[3]])]&/@Transpose[{Avec,Bmat,amp}])/PsiPsi];
		If[ResModes,amp*=Transpose[rescalelist];err2s*=Transpose[rescalelist^2]],
		(*below is for Greedy Algorithm of coefficients*)
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ BmatU;
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,\[Alpha]mat}];
		FixedCount=1;
		UnknownCount=1;
		ampFixed=Transpose[FixedCoef];
		ampUnknown=Transpose[amp];
		(* Reorder the fixed amplitudes and unknown amplitudes to make them have the initial order*)
		While[FixedCount<=Length[idxFixed]&&UnknownCount<=Length[idxUnknown],
		If[idxFixed[[FixedCount]]<=idxUnknown[[UnknownCount]], 
			AppendTo[ampTotal,ampFixed[[FixedCount]]];FixedCount++,
			AppendTo[ampTotal,ampUnknown[[UnknownCount]]];UnknownCount++
		]
		];
		While[FixedCount<=Length[idxFixed],AppendTo[ampTotal,ampFixed[[FixedCount]]];FixedCount++];
		While[UnknownCount<=Length[idxUnknown],AppendTo[ampTotal,ampUnknown[[UnknownCount]]];UnknownCount++];
		ampTotal=Transpose[ampTotal];
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		(*\[Rho]=Sqrt[Abs[((Conjugate[#[[1]]] . #[[2]])&/@Transpose[{\[Alpha]mat,amp}])/PsiPsiU]];*)
		\[Rho]=Sqrt[Abs[(((Conjugate[#[[2]]] . #[[1]])^2/(Conjugate[#[[2]]] . (#[[3]] . #[[2]])))&/@Transpose[{\[Alpha]mat,amp,BmatU}])/PsiPsiU]];
		\[Rho]T=Sqrt[Abs[(((Conjugate[#[[2]]] . #[[1]])^2/(Conjugate[#[[2]]] . (#[[3]] . #[[2]])))&/@Transpose[{Avec,ampTotal,Bmat}])/PsiPsi]];
		If[ResModes,amp*=Transpose[rescalelistUnknown];err2s*=Transpose[rescalelistUnknown^2]],
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ Bmat;
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,Avec}];	
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		\[Rho]=Sqrt[(Abs[Conjugate[#[[1]]] . (#[[3]])]^2/Abs[Conjugate[#[[3]]] . #[[2]] . (#[[3]])]&/@Transpose[{Avec,Bmat,amp}])/PsiPsi];
		If[ResModes,amp*=Transpose[rescalelist];err2s*=Transpose[rescalelist^2]]
		],
	(* Least Squares - DesignMatrix method *) 
		DesignMatrix,
		{indf,nsims,dm,b,IndTime,\[Omega]Info}=KRFDesignMatrix[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->OptionValue[NLmodesList]];
		rescalelist={};
		(*Note that the structure of the rescalelist in DesignMatrix is the transposed rescalelist of other method, e.g., Eigenvalue method and NormalEquation*)
		Switch[
		Head[TimeStride],
		Symbol,
		appends=Last@Reap[For[i=1,i<=indf,++i,
			dmCurr=dm;
			If[ResModes,
				AppendTo[rescalelist,Flatten[Exp[-Im[\[Omega]Info]*tr[[i]]]]]; 
				dmCurr*=rescalelist[[-1]];
			];		
			{u,w,v}=SingularValueDecomposition[SetPrecision[Transpose[dmCurr],prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[i]]]^2/Abs[Conjugate[acol] . (Bmat[[i]] . acol)]))/PsiPsi[[i]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			dm=Drop[dm,None,nsims]; b=Drop[b,nsims] (* remove data from latest time *)
		]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[ResModes,amp*=rescalelist;err2s*=rescalelist^2];
		If[returnsv,w:=appends[[4]]],
		Integer(* choose to fit with specific integer time stride *),
		IndCount=1;
		appends=Last@Reap[For[i=1,i<=indf,i+=TimeStride,
			dmCurr=dm;
			If[ResModes,
				AppendTo[rescalelist,Flatten[Exp[-Im[\[Omega]Info]*tr[[i]]]]]; 
				dmCurr*=rescalelist[[-1]];
			];
			{u,w,v}=SingularValueDecomposition[SetPrecision[Transpose[dmCurr],prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[IndCount]]]^2/Abs[Conjugate[acol] . (Bmat[[IndCount]] . acol)]))/PsiPsi[[IndCount]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			dm=Drop[dm,None,nsims*TimeStride]; b=Drop[b,nsims*TimeStride];(* remove data from latest time *)
			++IndCount 
		]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[ResModes,amp*=rescalelist;err2s*=rescalelist^2];
		If[returnsv,w:=appends[[4]]],
		List(* choose to fit with a time list *),
		dm=Drop[dm,None,nsims*(IndTime[[1]]-1)]; b=Drop[b,nsims*(IndTime[[1]]-1)];
		IndTimeDif=Differences[IndTime];
		IndCount=1;
		appends=Last@Reap[Do[
			dmCurr=dm;
			If[ResModes,
				AppendTo[rescalelist,Flatten[Exp[-Im[\[Omega]Info]*tr[[i]]]]]; 
				dmCurr*=rescalelist[[-1]];
			];
			{u,w,v}=SingularValueDecomposition[SetPrecision[Transpose[dmCurr],prec],Tolerance->Sqrt[tol]];
			winv=Map[If[#==0,0,1/#]&,(Diagonal[w])];
			acol=v . ((b . Conjugate[Take[u,All,Length[winv]]])winv);
			Sow[acol,amp];
			Sow[((Abs[v]^2) . (winv^2)),err2s];
			Sow[Sqrt[((Abs[Conjugate[acol] . Avec[[IndCount]]]^2/Abs[Conjugate[acol] . (Bmat[[IndCount]] . acol)]))/PsiPsi[[IndCount]]],\[Rho]];
			If[returnsv,Sow[w^2,w2]];
			If[IndCount<=Length[IndTimeDif],
			dm=Drop[dm,None,nsims*IndTimeDif[[IndCount]]]; b=Drop[b,nsims*IndTimeDif[[IndCount]]]]; (* remove data from latest time *)
			++IndCount
		,{i,IndTime}]];
		amp:=appends[[1]];
		err2s:=appends[[2]];
		\[Rho]:=appends[[3]];
		If[ResModes,amp*=rescalelist;err2s*=rescalelist^2];
		If[returnsv,w:=appends[[4]]]		
		],
	(* Least Squares - NormalEquation method *)
		NormalEquation,
		{indf,nsims,dm,b,IndTime,\[Omega]Info}=KRFDesignMatrix[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->OptionValue[NLmodesList]];
		rescalelist={};
		nk=Dimensions[dm][[1]];
		For[i=1,i<=nk,++i,
			conji=Conjugate[dm[[i]]];
			dots=Reverse[Take[Accumulate[Reverse[conji b],Method->"CompensatedSummation"],{-(indf-1)*nsims-1,-1,nsims}]];
			If[ResModes,AppendTo[rescalelist,Exp[-Im[\[Omega]Info[[i]]]*tr]]; dots*=rescalelist[[-1]]];
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
				If[ResModes, dots *=  Exp[-Im[\[Omega]Info[[i]]+\[Omega]Info[[j]]]*tr]];
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
		(*If we apply rescaling, select the rescale list basing on the FitTimeStride Options *)
		If[ResModes,
		Switch[Head[TimeStride],
			Symbol,
			Null,
			Integer,
			rescalelist=Transpose[(Transpose[rescalelist])[[#]]&/@Range[1,indf,TimeStride]],
			List,
			rescalelist=Transpose[(Transpose[rescalelist])[[#]]&/@IndTime]
		];];
		svd=SingularValueDecomposition[#,Tolerance->tol]& /@ SetPrecision[Normal[SparseArray[NMpos->Flatten[#]]]&/@Transpose[NMelem],prec];
		w=Flatten[Take[svd,All,{2}],1];v=Flatten[Take[svd,All,{3}],1];Clear[svd];
		winv=Map[If[#==0,0,1/#]&,(Diagonal[#]& /@ w),{2}];
		amp=#[[1]] . DiagonalMatrix[#[[2]]] . ConjugateTranspose[#[[1]]] . #[[3]]&/@Transpose[{v,winv,Transpose[Rvec]}];
		err2s=(Abs[#[[1]]]^2) . #[[2]]&/@Transpose[{v,winv}];
		\[Rho]=Sqrt[(Abs[Conjugate[#[[3]]] . #[[1]]]^2/Abs[Conjugate[#[[3]]] . #[[2]] . #[[3]]]&/@Transpose[{Avec,Bmat,amp}])/PsiPsi];
		If[ResModes,amp*=Transpose[rescalelist];err2s*=Transpose[rescalelist^2]]
	];
	If[Not[returnsv],
		{t,If[fullMismatch==True && Head[fixedModesGreedy]==List,\[Rho]T,\[Rho],\[Rho]],SetPrecision[amp,MachinePrecision],
			If[Head[NLlist]==List,
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,NLlist,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta],NLlist},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta],NLlist}],
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]}],
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]}]
				],
				SetPrecision[err2s,MachinePrecision],count},
		{t,If[fullMismatch==True && Head[fixedModesGreedy]==List,\[Rho]T,\[Rho],\[Rho]],SetPrecision[amp,MachinePrecision],
			If[Head[NLlist]==List,
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,NLlist,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta],NLlist},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta],NLlist}],
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]}],
				If[Head[fixedModesGreedy]==List,
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,{FixedModesG,paramGreedy}},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]},
					{{massratio,a,\[Theta],\[Phi]},SimModes,qnmp\[Delta],qnmm\[Delta]}]
				],
				SetPrecision[err2s,MachinePrecision],count,Diagonal[#]&/@w}
	]
]


Options[RestrictOverlap]=Union[{OmitModes->{{},{},{}}},Options[OverlapFit]];
(*This version of RestrictOverlap is compatible with fit obtained with greedy algorithm and nonlinear modes*)
RestrictOverlap[fit_List,opts:OptionsPattern[]]:=
Module[{prec=OptionValue[SVDWorkingPrecision],
	t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,omitp,omitm,omitnl,omitposp,omitposm,omitposnl,omitpos,
	Avec,Bmat,PsiPsi,count,qnmp,qnmm,err2,rescalelist,ResModes=OptionValue[RescaleModes],fittimeT0=OptionValue[T0],
	fittimeTFinal=OptionValue[TFinal],ind0,indf,tIndList,tRestrict,ampInd,ampRestrict,fitInfo,singularValue,fullMismatch=OptionValue[FullMismatch],
    qnmpSet, qnmmSet, fixedGreedyIndex, FixedModesG, QQNModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,
    nplus, nminus, idxUnknown, qnmpInd={},qnmmInd={}, NLInd, idxFixed={}, nUnknown, nFixed,
    AvecK, AvecU, Btrans, BmatK,BKelem={},BKpos={},BmatU,BUelem={},BUpos={},BmatC,BCelem={},BCpos={},
    BmatCconj, FixedCoef, rescalelistUnknown, rescalelistFixed, \[Alpha]mat, PsiPsiU,
    FixedCount, UnknownCount, ampFixed, ampUnknown, ampTotal={}, NLlist},
	{{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{qnmp,qnmm,QQNModes},{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}} = FitInfoStruct[fit];
	tIndList=TimeIndex[#]&/@t;
	If[fittimeT0<tIndList[[1]]||fittimeTFinal>tIndList[[-1]],Message[RestrictOverlap::TimeError,t[[1]],t[[-1]],KRFtime[[fittimeT0]],KRFtime[[fittimeTFinal]]];Abort[]];
	{omitp,omitm,omitnl}=OptionValue[OmitModes];
	(* Set QNM mode lists.  First Sort them. *)
(*	qnmp=SortBy[QNModesp,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];
	qnmm=SortBy[QNModesm,{Abs[#[[2]]]&,-#[[2]]&,#[[1]]&,#[[3]]&}];*)
	omitnl=DeleteDuplicates[Sort[#]&/@omitnl];
	omitposp=Table[If[MemberQ[Flatten[Position[qnmp,#]&/@omitp],i],0,1],{i,Length[qnmp]}];
	omitposm=Table[If[MemberQ[Flatten[Position[qnmm,#]&/@omitm],i],0,1],{i,Length[qnmm]}];
	omitposnl=Table[If[MemberQ[Flatten[Position[QQNModes,#]&/@omitnl],i],0,1],{i,Length[QQNModes]}];
	omitpos=Join[omitposp,omitposm,omitposnl];
	SetModeData[a,qnmpSet,qnmmSet,Evaluate@FilterRules[{opts},Options@SetModeData]];
	nplus=Length[qnmp];
	nminus=Length[qnmm];
	If[QQNModes=={},QQNModes=False];
	If[Length[FixedModesG]!=0,
	{Avec,Bmat,PsiPsi,tRestrict,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,fixedGreedyIndex,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->QQNModes],
	{Avec,Bmat,PsiPsi,tRestrict,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->QQNModes],
	{Avec,Bmat,PsiPsi,tRestrict,count,rescalelist}=ComputeInnerProducts[{massratio,a,\[Theta],\[Phi]},SimModes,qnmp,qnmm,
		T0->OptionValue[T0],TFinal->OptionValue[TFinal],TEnd->OptionValue[TEnd],RestrictToSimulationSubspace->OptionValue[RestrictToSimulationSubspace]
		,RescaleModes->OptionValue[RescaleModes],FitTimeStride->OptionValue[FitTimeStride],NLmodesList->QQNModes]	
	];
	Avec=SetPrecision[Avec,prec];Bmat=SetPrecision[Bmat,prec];PsiPsi=SetPrecision[PsiPsi,prec];
	If[Complement[tRestrict,t]!={},Message[RestrictOverlap::TimeStrideError,Complement[tRestrict,t]];Abort[]];
	ampInd=Flatten[Position[t,#]&/@tRestrict];
	ampRestrict=amp[[#]]&/@ampInd;	
	NLlist = QQNModes;
	(*If fit is obtained with greedy algorithm applied*)
	If[Length[FixedModesG]!=0,
	If[Head[NLlist]==List(*If there are nonlinear modes in the list*), 
		(*obtain the indices of QNMs and nonlinear modes *)
		idxUnknown=Range[1,nplus+nminus+Length[NLlist]];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		NLInd={#}&/@fixedGreedyIndex[[3]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus,fixedGreedyIndex[[3]]+nplus+nminus];
		idxFixed={#}&/@idxFixed;
		NLlist = Delete[NLlist,NLInd],
		(*If no nonlinear modes in the list*)
		idxUnknown=Range[1,nplus+nminus];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus];
		idxFixed={#}&/@idxFixed,
		idxUnknown=Range[1,nplus+nminus];
		qnmpInd={#}&/@fixedGreedyIndex[[1]];
		qnmmInd={#}&/@fixedGreedyIndex[[2]];
		idxFixed=Join[fixedGreedyIndex[[1]],fixedGreedyIndex[[2]]+nplus];
		idxFixed={#}&/@idxFixed
	];
	qnmp = Delete[qnmp,qnmpInd];
	qnmm = Delete[qnmm,qnmmInd];
	idxUnknown=Delete[idxUnknown,idxFixed];
	idxFixed=Flatten[idxFixed];
	nUnknown=Length[idxUnknown];
	nFixed=Length[idxFixed];	
	(*Seperate unknown and known parts from the A vector *)
	AvecK=Transpose[Part[Transpose[Avec],idxFixed]];
	AvecU=Transpose[Part[Transpose[Avec],idxUnknown]];
	(*Seperate B fixed from B matrix*)
	Btrans=Transpose[Bmat,{3,1,2}];
	Do[AppendTo[BKelem,Btrans[[idxFixed[[i]],idxFixed[[j]]]]];AppendTo[BKpos,{i,j}];
		If[j>i,AppendTo[BKelem,Conjugate[Btrans[[idxFixed[[i]],idxFixed[[j]]]]]];AppendTo[BKpos,{j,i}]]
		,{i,nFixed},{j,i,nFixed}];
	BmatK=Normal[SparseArray[BKpos->Flatten[#]]]&/@Transpose[BKelem];
	(*Seperate B unknown from B matrix *)
	Do[AppendTo[BUelem,Btrans[[idxUnknown[[i]],idxUnknown[[j]]]]];AppendTo[BUpos,{i,j}];
		If[j>i,AppendTo[BUelem,Conjugate[Btrans[[idxUnknown[[i]],idxUnknown[[j]]]]]];AppendTo[BUpos,{j,i}]]
		,{i,nUnknown},{j,i,nUnknown}];
	BmatU=Normal[SparseArray[BUpos->Flatten[#]]]&/@Transpose[BUelem];	
	(*Seperate B fixed and unknown overlap part from B matrix *)
	Do[AppendTo[BCelem,Btrans[[idxFixed[[i]],idxUnknown[[j]]]]];AppendTo[BCpos,{i,j}]
		,{i,nFixed},{j,nUnknown}];
	BmatC=Normal[SparseArray[BCpos->Flatten[#]]]&/@Transpose[BCelem];
	BmatCconj=ConjugateTranspose[#]&/@BmatC;
	(*Set up fixed Coefficients. *)
	FixedCoef=Transpose[(Table[1,{i,Length[Avec]}]*#[[2]])&/@FixedModesG];
	(* Get the fixed and unknown parts of the rescalelist. If rescaling is applied, the 
	fixed coefficients should also be changed to match the rescaling *)
	If[ResModes,
		rescalelistUnknown=Part[rescalelist,idxUnknown];
		rescalelistFixed=Part[rescalelist,idxFixed];
		FixedCoef/=Transpose[rescalelistFixed];
		ampRestrict/=Transpose[rescalelistUnknown]];
	(*Calculate the unknown part of PsiPsi*)
	\[Alpha]mat=(#[[1]]-#[[2]] . #[[3]])&/@Transpose[{AvecU,BmatCconj,FixedCoef}];
	PsiPsiU=(#[[1]]-2*Re[Conjugate[#[[4]]] . #[[3]]]+Conjugate[#[[4]]] . (#[[2]] . #[[4]]))&/@Transpose[{PsiPsi,BmatK,AvecK,FixedCoef}];

	FixedCount=1;
	UnknownCount=1;
	ampFixed=Transpose[FixedCoef];
	ampUnknown=Transpose[ampRestrict];
	(* Reorder the fixed amplitudes and unknown amplitudes to make them have the initial order*)
	While[FixedCount<=Length[idxFixed]&&UnknownCount<=Length[idxUnknown],
		If[idxFixed[[FixedCount]]<=idxUnknown[[UnknownCount]], 
			AppendTo[ampTotal,ampFixed[[FixedCount]]];FixedCount++,
			AppendTo[ampTotal,ampUnknown[[UnknownCount]]];UnknownCount++
		]
	];
	While[FixedCount<=Length[idxFixed],AppendTo[ampTotal,ampFixed[[FixedCount]]];FixedCount++];
	While[UnknownCount<=Length[idxUnknown],AppendTo[ampTotal,ampUnknown[[UnknownCount]]];UnknownCount++];
	ampTotal=Transpose[ampTotal];
	If[fullMismatch==True,
		\[Rho]=Sqrt[Abs[(((Conjugate[(#[[2]]omitpos)] . #[[1]])^2/(Conjugate[(#[[2]]omitpos)] . (#[[3]] . (#[[2]]omitpos))))&/@Transpose[{Avec,ampTotal,Bmat}])/PsiPsi]],
		\[Rho]=Sqrt[Abs[(((Conjugate[(#[[2]]omitpos)] . #[[1]])^2/(Conjugate[(#[[2]]omitpos)] . (#[[3]] . (#[[2]]omitpos))))&/@Transpose[{\[Alpha]mat,ampRestrict,BmatU}])/PsiPsiU]]			
	];
	];
	(*If fit is not obtained with greedy algorithm applied*)
	If[Length[FixedModesG]==0,
		If[ResModes,ampRestrict/=Transpose[rescalelist]];
		\[Rho]=Sqrt[(Abs[Conjugate[#[[1]]] . (#[[3]]omitpos)]^2/Abs[Conjugate[#[[3]]omitpos] . #[[2]] . (#[[3]]omitpos)]&/@Transpose[{Avec,Bmat,ampRestrict}])/PsiPsi]
	];
	{tRestrict,SetPrecision[\[Rho],MachinePrecision]}
]


Options[FitAmplitudesTable]=Options[SetModeData];
FitAmplitudesTable[fit_List,fittime_?NumberQ,opts:OptionsPattern[]]:=
Module[{t,amp,ind,ind2,QNModesp,QNModesm,qnmfixedp,qnmfixedm,massratio,a,\[Theta],\[Phi],SimModes,\[Rho],QNMamp,count,
        l,m,n,lp,mp,np,qnmp,qnmm,qqnm,mode,err2s,cerr=0,
        qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG,QQNModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,err2,
        j,massratioFixed,QNModespFixed,QNModesmFixed,QQNModesFixed={},NLmodesExist,QQNModes\[Omega]Unfixed,QQNModes\[Omega]Fixed,
        \[Omega]NL,k,fixedAmp,qqnmIndex},
    {{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{QNModesp,QNModesm,QQNModes},{t,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}}= FitInfoStruct[fit];
	If[QQNModes=={},NLmodesExist=False,NLmodesExist=True,NLmodesExist=True];
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
	ind=If[fittime>=t[[-1]],Length[t],If[fittime<=t[[1]],1,SequencePosition[t,{x_/;x>=fittime},1][[1,1]]]];
	ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
	SetModeData[a,qnmpSet,qnmmSet,Evaluate@FilterRules[{opts},Options@SetModeData]];
	Do[{l,m}=Clm;
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
		cerr+=Conjugate[mode] . mode,
	{Clm,SimModes}];
	err2s=Sqrt[err2[[ind]]*Abs[cerr]/(2(Length[SimModes]count[[ind]]-Length[QNModespUnfixed]-Length[QNModesmUnfixed]-Length[QQNModesUnfixed]))];
	TextGrid[Join[
		{{"Mode","Amplitude","Phase/\[Pi]","\[Sigma](Amp)","\[Sigma](Phase)/\[Pi]"}},
		Flatten[#,{1,2}]&/@Transpose[{Join[{ToString[#]<>"+"}&/@QNModespUnfixed,{ToString[#]<>"-"}&/@QNModesmUnfixed,{ToString[#]}&/@QQNModesUnfixed],
										AbsArg[Chop[#]]/{1,\[Pi]}&/@amp[[ind]],
										Transpose[{err2s,If[Im[#]!=0,1,#/\[Pi]]&/@ArcSin[err2s/Abs[amp[[ind]]]]}]}]
			],Spacings->2]
]




Options[RelativeAmplitudes]=Options[SetModeData];
RelativeAmplitudes[fit_List,lm_List,fittime_,opts:OptionsPattern[]]:=
Module[{t,\[Rho],amp,massratio,a,\[Theta],\[Phi],SimModes,QNModesp,QNModesm,count,ind,l,m,qnmp,qnmm,qqnm,lp,mp,np,j,relamps={},
	qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG,QQNModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,err2,
	NLmodesExist,massratioFixed,QNModespFixed,QNModesmFixed,QQNModesFixed={},QQNModes\[Omega]Unfixed,QQNModes\[Omega]Fixed,\[Omega]NL},
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
	ind=If[NumberQ[fittime],
		If[fittime>=t[[-1]],Length[t],If[fittime<=t[[1]],1,SequencePosition[t,{x_/;x>=fittime},1][[1,1]]]],
		Switch[Head[fittime],
			Symbol,If[fittime==All,Range[1,Length[t]],Message[RelativeAmplitudes::Abort,fittime];Abort[],Message[RelativeAmplitudes::Abort,fittime];Abort[]],
			_,Message[RelativeAmplitudes::Abort,fittime];Abort[]
			]
		];
	{l,m}=lm;
	qnmp=SpheroidalHarmonicModes[lm,QNModespUnfixed];
	qnmm=SpheroidalHarmonicModes[lm,QNModesmUnfixed];
	If[UseSpheroidalExpansion,
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,Qlmnp]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] t/massratio]]],
		{Qlmnp,qnmp}];
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,Qlmnp]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] t/massratio]]],
		{Qlmnp,qnmm}];
		If[NLmodesExist, (*if there are quadratic QNMs*)
				qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Unfixed];
				qqnmIndex=Flatten[Position[QQNModes\[Omega]Unfixed,#]&/@qqnm];
				For[j=1,j<=Length[qqnm],j++,
					{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
					AppendTo[relamps, Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*amp[[ind,FindNonlinearIndex[QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,QQNModesUnfixed[[qqnmIndex[[j]]]]]]]*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* t/massratio]]]	
				];
		]
		,
	(* Use this version to ignore Spheroidal Expansion expansion coefficients *)
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[KroneckerDelta[l,lp]amp[[ind,QNMpIndex[QNModespUnfixed,QNModesmUnfixed,Qlmnp]]]Exp[-I KRF\[Omega][lp,mp,np] t/massratio]]],
		{Qlmnp,qnmp}];
		Do[{lp,mp,np}=Qlmnp;
			AppendTo[relamps,Abs[KroneckerDelta[l,lp]amp[[ind,QNMmIndex[QNModespUnfixed,QNModesmUnfixed,Qlmnp]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] t/massratio]]],
		{Qlmnp,qnmm}];
		If[NLmodesExist, (*if there are quadratic QNMs*)
				qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Unfixed];
				qqnmIndex=Flatten[Position[QQNModes\[Omega]Unfixed,#]&/@qqnm];
				For[j=1,j<=Length[qqnm],j++,
					{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
					AppendTo[relamps, Abs[WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*amp[[ind,FindNonlinearIndex[QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,QQNModesUnfixed[[qqnmIndex[[j]]]]]]]*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* t/massratio]]]	
				];
		]
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
	(*Comparing to the version 1.0.0. It will return the information of nonlinear modes and fixed modes*)
	{fits[[1,1,1,1]],mm,amp,Join[{bhp},Flatten[Take[fits,1,1,1,{4},{2,-1}],4]],err,cnt}
]


Options[MaximizeOverlap]=Union[{},Options[OverlapFit],Options[FindMaximum]];
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
			FindMaximum[
			(*OverlapFit[\[Delta]f,\[Chi]f,0,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],
						{{\[Delta]f,\[Delta]g,1.00001\[Delta]g},{\[Chi]f,\[Chi]g,1.00001\[Chi]g}},*)						
						{OverlapFit[\[Delta]f,\[Chi]f,0,0,sm,qnmp,qnmm,T0->ti,TFinal->ti,ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]][[2,1]],
						\[Delta]g-0.02<\[Delta]f<\[Delta]g+0.02,\[Chi]g-0.02<\[Chi]f<\[Chi]g+0.02},
						(*The restrictions \[Delta]g-0.02<\[Delta]f<\[Delta]g+0.02 and \[Chi]g-0.02<\[Chi]f<\[Chi]g+0.02 mean to help prevent the nonlinear fitting routine
						from finding unrealistic solution of {\[Delta]f,\[Chi]f}. It is useful during the nonlinear greedy approach.
						Because that kinds of fitting tends to have more noise in the parameter space. *)
						{{\[Delta]f,\[Delta]g,1.00001*\[Delta]g},{\[Chi]f,\[Chi]g,1.00001*\[Chi]g}},
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


Options[RemnantParameterSpaceMaxOverlap]=Union[{FitAngle->False},Options[OverlapFit],Options[FindMaximum]];
RemnantParameterSpaceMaxOverlap[rps_List,opts:OptionsPattern[]]:=
Module[{Nt,params,i,pos,minparam,fit\[Theta]=OptionValue[FitAngle],maxo,ret={},
fixedModesGreedy=OptionValue[FixedModesGreedy],qnmp,qnmm,FixedModes,FixedQNM={},tRps=rps[[1]],
tIndList,fittimeT0=OptionValue[T0],fittimeTFinal=OptionValue[TFinal],t,tLength,
TimeStride=OptionValue[FitTimeStride],TimePos,IndTime,IndTimeRps,tInd},
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
	qnmp=rps[[4,3]];
	qnmm=rps[[4,4]];
	If[Head[fixedModesGreedy]==List,
		FixedModes=#[[1]]&/@fixedModesGreedy[[1]];
		If[Length[#]==4,AppendTo[FixedQNM,#]]&/@FixedModes;
		If[#[[4]]==1,AppendTo[qnmp,#[[1;;3]]]
		,AppendTo[qnmm,#[[1;;3]]]]&/@FixedQNM;
		qnmp=DeleteDuplicates[qnmp];
		qnmm=DeleteDuplicates[qnmm]
	];
	Nt=Length[t];
	For[i=1,i<=Nt,++i,
		tInd=IndTimeRps[[i]];
		pos=Flatten[Position[rps[[2,tInd]],Min[rps[[2,tInd]]]],1];
		(*minparam=Flatten[rps[[4,1,pos[[1]],pos[[2]],pos[[3]]]]];*)
		minparam=rps[[4,1,pos[[1]],pos[[2]],pos[[3]]]];
		maxo=If[fit\[Theta],
				MaximizeOverlap[minparam[[1]],minparam[[2]],minparam[[3]],
								TimeIndex[t[[i]]],rps[[4,2]],qnmp,qnmm,
								Evaluate@FilterRules[{opts},Options@OverlapFit]],			
				MaximizeOverlap[minparam[[1]],minparam[[2]],
								TimeIndex[t[[i]]],rps[[4,2]],qnmp,qnmm,
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


Options[MaxOverlapSequenceAmplitudes]=Union[{},Options[OverlapFit]];
MaxOverlapSequenceAmplitudes[mos_List,SimModesInput_List,QNModespInput_List,QNModesmInput_List,opts:OptionsPattern[]]:=
Module[{fittime,massratio,a,\[Theta],ind,ind2,fit,amp,err2,count,
		i,l,m,mode,qnmp,qnmm,j,lp,mp,np,cerr,t={},amps={},err2s={},
		qnmpSet,qnmmSet,fixedGreedyIndex={},FixedModesG={},QQNModes,\[Rho],\[Phi],QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed={},
		massratioFixed,QNModespFixed,QNModesmFixed,QQNModesFixed={},QQNModes\[Omega]Unfixed,QQNModes\[Omega]Fixed,NLmodesExist,
		qqnm,qqnmIndex,\[Omega]NL,k,fixedAmp,QNModesp,QNModesm,SimModes,tTemp},
	For[ind=1,ind<=Length[mos],++ind,
		fittime=mos[[ind,1]];
		ind2=If[fittime>=KRFtime[[-1]],Length[KRFtime],If[fittime<=KRFtime[[1]],1,SequencePosition[KRFtime,{x_/;x>=fittime},1][[1,1]]]];
		{massratio,a,\[Theta]}=mos[[ind,{3,4,5}]];		
		fit=OverlapFit[massratio,a,\[Theta],0,DeleteDuplicates[SimModesInput],
			DeleteDuplicates[QNModespInput],DeleteDuplicates[QNModesmInput],
			T0->ind2,TFinal->ind2,
			ReturnSingularValues->False,Evaluate@FilterRules[{opts},Options@OverlapFit]];
		{{qnmpSet,qnmmSet,fixedGreedyIndex,FixedModesG},{QNModesp,QNModesm,QQNModes},{tTemp,\[Rho],amp,{{massratio,a,\[Theta],\[Phi]},SimModes,QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed},err2,count}} = FitInfoStruct[fit];
		amp=amp[[1]];
		err2=err2[[1]];
		count=count[[1]];
		If[QQNModes=={},NLmodesExist=False,NLmodesExist=True,NLmodesExist=True];
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
		cerr=0;
		For[i=1,i<=Length[SimModes],++i,
			{l,m}=SimModes[[i]];
			mode=Table[0,Length[KRFtime]];
			(*contribution to waveform from unfixed QNM plus*)
			qnmp=SpheroidalHarmonicModes[{l,m},QNModespUnfixed];
			For[j=1,j<=Length[qnmp],++j,
				{lp,mp,np}=qnmp[[j]];
				mode+=If[UseSpheroidalExpansion,
					WignerD[{l,-m,-mp},0,\[Theta],0]amp[[QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]KRFYS[l,lp,mp,np]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio],
        			(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
        			KroneckerDelta[l,lp]amp[[QNMpIndex[QNModespUnfixed,QNModesmUnfixed,qnmp[[j]]]]]Exp[-I KRF\[Omega][lp,mp,np] KRFtime/massratio]
        		];
			];
			(*contribution to waveform from unfixed QNM minus*)
			qnmm=SpheroidalHarmonicModes[{l,m},QNModesmUnfixed];
			For[j=1,j<=Length[qnmm],++j,
				{lp,mp,np}=qnmm[[j]];
				mode+=If[UseSpheroidalExpansion,
					(-1)^(l+lp)WignerD[{l,-m,-mp},\[Phi],\[Theta],0]amp[[QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Conjugate[KRFYS[l,lp,-mp,np]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio],
					(* Use this version of mode+= to ignore the Spheroidal Harmonic expansion coefficients *)
					KroneckerDelta[l,lp]amp[[QNMmIndex[QNModespUnfixed,QNModesmUnfixed,qnmm[[j]]]]]Exp[I Conjugate[KRF\[Omega][lp,-mp,np]] KRFtime/massratio]
				];
			];
			If[NLmodesExist, (*if there are quadratic QNMs*)
				qqnm=SpheroidalHarmonicModesNL[{l,m},QQNModes\[Omega]Unfixed];
				qqnmIndex=Flatten[Position[QQNModes\[Omega]Unfixed,#]&/@qqnm];
				For[j=1,j<=Length[qqnm],j++,
					{lp,mp,\[Omega]NL}=qqnm[[j]];(* nonlinear modes*)
					mode+=WignerD[{l,-m,-mp},\[Phi],\[Theta],0]*amp[[FindNonlinearIndex[QNModespUnfixed,QNModesmUnfixed,QQNModesUnfixed,QQNModesUnfixed[[qqnmIndex[[j]]]]]]]*KroneckerDelta[l,lp]Exp[-I \[Omega]NL* KRFtime/massratio]	
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
        	mode=Take[Drop[mode,ind2-1],count];
        	cerr+=Conjugate[mode] . mode;
        ];
        err2=Sqrt[err2*Abs[cerr]/(2(Length[SimModes]count-Length[QNModespUnfixed]-Length[QNModesmUnfixed]-Length[QQNModesUnfixed]))];        
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
