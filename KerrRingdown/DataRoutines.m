(* ::Package:: *)

(* ::Chapter:: *)
(*Data Routines Package*)


(* ::Section::Closed:: *)
(*Begin DataRoutines Package*)


BeginPackage["DataRoutines`"]


(* ::Section::Closed:: *)
(*Documentation of External Functions*)


(* ::Subsection::Closed:: *)
(*SXS Data Routines*)


SXSWaveform::usage=
"SXSWaveform[\!\(\*
StyleBox[\"sxsdir\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"Next\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)]\n" <>
"\t \!\(\*
StyleBox[\"sxsdir\",\nFontSlant->\"Italic\"]\) : String containing the full path to the directory containing simulation data\n"<>
"\t \!\(\*
StyleBox[\"Next\",\nFontSlant->\"Italic\"]\) : Waveform extrapolation to use:\n"<>
"\t\t 0 : Outermost extraction radius.\n" <>
"\t\t 2 : Extrapolated N2.\n"<>
"\t\t 3 : Extrapolated N3.\n"<>
"\t\t 4 : Extrapolated N4.\n"<>
"\t \!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\) : Spin-weighted spherical harmonic indices of mode to be read into memeory.\n\n"<>
"Options:\n"<>
"\t WaveformType : Defaults to Metric.\n"<>
"\t\t Psi4 : Read from rMPsi4 files\n"<>
"\t\t Metric : Read from rhoverM files\n"<>
"\t FrameType : Defaults to Raw.\n"<>
"\t\t Raw : Read from {WaveformType}_Asymptotic_GeometricUnits.h5\n"<>
"\t\t CoM : Read from {WaveformType}_Asymptotic_GeometricUnits_CoM.h5"


SXSCCEWaveform::usage=
"SXSCCEWaveform[\!\(\*
StyleBox[\"sxsdir\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"RNext\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\)]\n" <>
"\t \!\(\*
StyleBox[\"sxsdir\",\nFontSlant->\"Italic\"]\) : String containing the full path to the directory containing simulation data\n"<>
"\t \!\(\*
StyleBox[\"RNext\",\nFontSlant->\"Italic\"]\) : CCE radii of a simulation or waveform extrapolation order to use.\n"<>
"\t \!\(\*
StyleBox[\"l\",\nFontSlant->\"Italic\"]\),\!\(\*
StyleBox[\"m\",\nFontSlant->\"Italic\"]\) : Spin-weighted spherical harmonic indices of mode to be read into memeory.\n\n"<>
"Options:\n"<>
"\t WaveformType : Defaults to Metric.\n"<>
"\t\t Psi4 : Read from rMPsi4 files\n"<>
"\t\t News : Read from r2News files\n"<>
"\t\t Metric : Read from rhoverM files\n"<>
"\t FrameType : Defaults to CoM.\n"<>
"\t\t CoM : Read from {WaveformType}..._CoM.h5\n"<>
"\t\t Mem : Read from {WaveformType}..._CoM_Mem.h5\n"<>
"\t Extrapolated : Defaults to False.\n"<>
"\t\t True : Extrapolated data is used. Read from {WaveformType}_Extrapolated_N{RNext}_CoM.h5\n"<>
"\t\t False : Extrapolated data is not used. CCE radii of a simulation is read in from RNext. Read from {WaveformType}_BondiCce_R{RNext}_CoM.h5\n"<>
"\t Superrest : Defaults to False.\n"<>
"\t\t True : Signal has been transformed to Super rest frame. Read from {WaveformType}..._CoM_Bondi.h5\n"<>
"\t\t False : Signal has not been tranformed to Super rest frame. Read from {Waveforminfo}..._CoM.h5\n"


SXSFinalProperties::usage=""


(* ::Subsection::Closed:: *)
(*Reserved Globals*)


Protect[WaveformType,Psi4,Metric,News,FrameType,Raw,CoM,ReM,Extrapolated,Superrest,Mem];


Begin["`Private`"]


(* ::Section::Closed:: *)
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


(* ::Section::Closed:: *)
(*BHPT Data Routines*)


Options[BHPTWaveform] = {WaveformType->Psi4,DataRange->All};
BHPTWaveform[h5name_String,l_Integer,m_Integer,OptionsPattern[]]:= 
Module[{mname,lname,Yname,WFname,rawdat,
        range=OptionValue[DataRange],
        wftype=SymbolName[OptionValue[WaveformType]]},
   mname = If[m<0,"-"<>ToString[Abs[m]],ToString[m]];
   lname = ToString[l];
   Yname = "Y_l"<>lname<>"_m"<>mname<>".dat";
   WFname = Switch[wftype,"Psi4","/DmuPsi4/",
                          "Metric","/DhOvermu/",
                          _,Abort[]];
   Off[Import::dataset];
   rawdat=Import[h5name,{"HDF5","Datasets",{WFname<>Yname}}];
   On[Import::dataset];
   If[rawdat==$Failed,Return[$Failed]];
   rawdat=Take[rawdat,range];
   {Flatten[Take[rawdat,All,{1}]],Function[x,x[[1]]+I x[[2]]]/@Drop[rawdat,0,1]}
]


(* ::Section::Closed:: *)
(*End of DataRoutines Package*)


End[] (* `Private` *)


EndPackage[]
