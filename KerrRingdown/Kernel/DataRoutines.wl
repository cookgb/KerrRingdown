(* ::Package:: *)

(* ::Chapter:: *)
(*Data Routines Package*)


(* ::Section::Closed:: *)
(*Begin DataRoutines Package*)


BeginPackage["KerrRingdown`"]


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*SXS Data Routines*)


SXSWaveform::usage=
"SXSWaveform[sxsdir,Next,l,m] "<>
"Read in waveform saved in spherical harmonic modes {l, m} "<>
"from the SXS catalogue source file in directory sxsdir."


SXSCCEWaveform::usage=
"SXSCCEWaveform[sxsdir,RNext,l,m] "<>
"Read in waveform saved in spherical harmonic modes {l, m} "<>
"from the Ext-CCE Waveform source file in directory sxsdir."


SXSFinalProperties::usage=
"SXSFinalProperties[\[Delta],\!\(\*SubscriptBox[\(\[Chi]\), \(x\)]\),\!\(\*SubscriptBox[\(\[Chi]\), \(y\)]\),\!\(\*SubscriptBox[\(\[Chi]\), \(z\)]\),\!\(\*SubscriptBox[\(v\), \(x\)]\),\!\(\*SubscriptBox[\(v\), \(y\)]\),\!\(\*SubscriptBox[\(v\), \(z\)]\)] "<>
"Gives the remnant Black Hole properties {\[Delta], \[Chi], \[Theta], \[Phi]}, which "<>
"can be used for further fitting, based on the information of "<>
"remnant parameters provided by the waveform."


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


(*Protect[WaveformType,Psi4,Metric,News,FrameType,Raw,CoM,ReM,Extrapolated,Superrest,Mem];*)


Protect[T0,TEnd,TFinal,Width,DataType,SXS,SXSCCE,SXSRNext,WaveformType,Psi4,Metric,News,
		FrameType,Extrapolated,Superrest,CoM,Mem,\[Delta]f,\[Chi]f,\[Theta]f,FitAngle,SVDWorkingPrecision,InitialGuess,
		CheckUpdate,Errors,UseLeastSquares,NormalEquation,ReturnSingularValues,
		Log10MisMatchRange,OmitModes,RestrictToSimulationSubspace,RotateFrame,FitTimeStride,RescaleModes];


Begin["`Private`"]


(* ::Section::Closed:: *)
(*SXS Data Routines *)


Options[SXSWaveform] = {WaveformType->Metric,DataRange->All,FrameType->CoM};
SXSWaveform[sxsdir_String,Next_Integer,l_Integer,m_Integer,OptionsPattern[]]:= 
Module[{mname,lname,Yname,Gname,WFname,h5name,rawdat,
        range=OptionValue[DataRange],
        wftype=SymbolName[OptionValue[WaveformType]],dtype=SymbolName[OptionValue[FrameType]]},
   mname = If[m<0,"-"<>ToString[Abs[m]],ToString[m]];
   lname = ToString[l];
   Yname = "Y_l"<>lname<>"_m"<>mname<>".dat";
   Gname = If[Next==0,"OutermostExtraction.dir/",
              If[Next==2,"Extrapolated_N2.dir/",
                 If[Next==3,"Extrapolated_N3.dir/",
                    If[Next==4,"Extrapolated_N4.dir/"]]]];
   WFname = Switch[wftype,"Psi4","/rMPsi4",
                          "Metric","/rhOverM",
                          _,Abort[]];
   h5name=Switch[dtype,"Raw",WFname<>"_Asymptotic_GeometricUnits.h5",
                       "CoM",WFname<>"_Asymptotic_GeometricUnits_CoM.h5",
                       _,Abort[]];
   Off[Import::dataset];
   rawdat=Import[sxsdir<>h5name,{"HDF5","Datasets",{Gname<>Yname}}];
   On[Import::dataset];
   If[rawdat==$Failed,Return[$Failed]];
   rawdat=Take[rawdat,range];
   {Flatten[Take[rawdat,All,{1}]],Function[x,x[[1]]+I x[[2]]]/@Drop[rawdat,0,1]}
]


Options[SXSCCEWaveform] = {WaveformType->Metric,DataRange->All,FrameType->CoM,Extrapolated->False,Superrest->False};
SXSCCEWaveform[sxsdir_String,RNext_Integer,l_Integer,m_Integer,OptionsPattern[]]:= 
Module[{mname,lname,Yname,Gname,WFname,h5name,rawdat,
        range=OptionValue[DataRange],
        wftype=SymbolName[OptionValue[WaveformType]],dtype=SymbolName[OptionValue[FrameType]]},
   mname = If[m<0,"-"<>ToString[Abs[m]],ToString[m]];
   lname = ToString[l];
   Yname = "Y_l"<>lname<>"_m"<>mname<>".dat";
   If[OptionValue[Extrapolated],
      Gname = "_Extrapolated_N"<>ToString[RNext],
      Gname = "_BondiCce_R"<>If[RNext>99,"0",If[RNext>9,"00","000"]]<>ToString[RNext]
      ];
   WFname = Switch[wftype,"Psi4","/rMPsi4",
                          "News","/r2News", 
                          "Metric","/rhOverM",
                          _,Abort[]];
   WFname=WFname<>Gname;    
   h5name=Switch[dtype,"CoM",WFname<>"_CoM",
                       "Mem",WFname<>"_CoM_Mem",
                       _,Abort[]];
   h5name=h5name<>If[OptionValue[Superrest],"_Bondi",""]<>".h5";                  
   Off[Import::dataset];
   rawdat=Import[sxsdir<>h5name,{"HDF5","Datasets",{Yname}}];
   On[Import::dataset];
   If[rawdat==$Failed,Return[$Failed]];
   rawdat=Take[rawdat,range];
   {Flatten[Take[rawdat,All,{1}]],Function[x,x[[1]]+I x[[2]]]/@Drop[rawdat,0,1]}
]


(* Assumes Eadm, J, and P are given as dimensionful quantities in units where
   the mass scale is 1 *)
(* Note: in metadata.txt, the coordinate velocity may not be P *)
(*SXSFinalProperties[Eadm_?NumberQ,jx_?NumberQ,jy_?NumberQ,jz_?NumberQ,
                   px:_?NumberQ:0,py:_?NumberQ:0,pz:_?NumberQ:0]:=
Module[{\[Gamma],j={jx,jy,jz},p={px,py,pz},a},
   \[Gamma]=Eadm/Sqrt[Eadm^2-p . p];
   a=\[Gamma]^2 j/Eadm^2-If[p . p==0,0,\[Gamma](\[Gamma]-1)p (j . p)/(p . p Eadm^2)];
   {Eadm/\[Gamma],
   Sqrt[a . a],
   ArcTan[a[[3]],Sqrt[a[[1]]^2+a[[2]]^2]],
   If[a[[1]]==a[[2]]==0,0.`,ArcTan[a[[1]],a[[2]]]]}
]*)

(*\[Delta] is Remnant mass. \[Chi] is dimensionless spin of BH remnant. v is coordinate velocity of remnant. 
It will give out remnant mass and dimensionless spin in rest frame of the BH. 
*)
SXSFinalProperties[\[Delta]_?NumberQ,\[Chi]x_?NumberQ,\[Chi]y_?NumberQ,\[Chi]z_?NumberQ,
                   vx:_?NumberQ:0,vy:_?NumberQ:0,vz:_?NumberQ:0]:=
Module[{\[Gamma],\[Gamma]\[Chi],\[Chi]={\[Chi]x,\[Chi]y,\[Chi]z},v={vx,vy,vz},\[Chi]p},
   \[Gamma]=1/Sqrt[1-v . v];
   \[Gamma]\[Chi]=1/Sqrt[1/2(1+Sqrt[1-\[Chi] . \[Chi]]+\[Gamma]^2(\[Chi] . \[Chi]-(v . \[Chi])^2)/(1+Sqrt[1-\[Chi] . \[Chi]]))];
   \[Chi]p=\[Gamma] \[Gamma]\[Chi]^2 (\[Chi]-\[Gamma] v (v . \[Chi])/(\[Gamma]+1));
   {\[Delta]/\[Gamma]\[Chi],
   Sqrt[\[Chi]p . \[Chi]p],
   ArcTan[\[Chi]p[[3]],Sqrt[\[Chi]p[[1]]^2+\[Chi]p[[2]]^2]],
   If[\[Chi]p[[1]]==\[Chi]p[[2]]==0,0.`,ArcTan[\[Chi]p[[1]],\[Chi]p[[2]]]]}
]


(* ::Section::Closed:: *)
(*End of DataRoutines Package*)


End[] (* `Private` *)


EndPackage[]
