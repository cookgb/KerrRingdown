(* ::Package:: *)

(* ::Chapter:: *)
(*SXS Data Routines Package*)


(* ::Section:: *)
(*Begin SXSDataRoutines Package*)


BeginPackage["SXSDataRoutines`"]


(* ::Section:: *)
(*Documentation of External Functions*)


(* ::Subsection:: *)
(*SXS Data Routines*)


(* ::Section:: *)
(*SXS Data Routines *)


Options[SXSWaveform] = {WaveformType->Metric,DataRange->All,FrameType->Raw};
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
                       "ReM",WFname<>"_Asymptotic_GeometricUnits_ReM.h5",
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

(*\[Delta] is Remnant mass. j is dimensionless spin of BH remnant. v is coordinate velocity of remnant. 
It will give out remnant mass and dimensionless spin in rest frame of the BH. 
*)
SXSFinalProperties[\[Delta]_?NumberQ,jx_?NumberQ,jy_?NumberQ,jz_?NumberQ,
                   vx:_?NumberQ:0,vy:_?NumberQ:0,vz:_?NumberQ:0]:=
Module[{\[Gamma],\[Gamma]\[Chi],j={jx,jy,jz},v={vx,vy,vz},a},
   \[Gamma]=1/Sqrt[1-v . v];
   \[Gamma]\[Chi]=1/Sqrt[1-(v . j)^2/(2 (1+Sqrt[1-j . j]))];
   a=\[Gamma]\[Chi]^2 (j-\[Gamma] v (v . j)/(\[Gamma]+1));
   {\[Delta]/\[Gamma]\[Chi],
   Sqrt[a . a],
   ArcTan[a[[3]],Sqrt[a[[1]]^2+a[[2]]^2]],
   If[a[[1]]==a[[2]]==0,0.`,ArcTan[a[[1]],a[[2]]]]}
]
