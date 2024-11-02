(* ::Package:: *)

BeginPackage["KerrRingdown`"];
(* Declare your packages public symbols here. *)


Begin["`Private`"];
(* Define your public and private symbols here. *)
If[Catch[Quiet[Check[Get["DataRoutinesLocal`"],Throw[True]],Get::noopen]],Get["DataRoutines`"],Print["Get Local DataRoutines"],Print["Get Local DataRoutines"]];
If[Catch[Quiet[Check[Get["ReadWaveformsLocal`"],Throw[True]],Get::noopen]],Get["ReadWaveforms`"],Print["Get Local ReadWaveforms"],Print["Get Local ReadWaveforms"]];
If[Catch[Quiet[Check[Get["QNMModeDataLocal`"],Throw[True]],Get::noopen]],Get["QNMModeData`"],Print["Get Local QNMModeData"],Print["Get Local QNMModeData"]];
If[Catch[Quiet[Check[Get["UtilityRoutinesLocal`"],Throw[True]],Get::noopen]],Get["UtilityRoutines`"],Print["Get Local UtilityRoutines"],Print["Get Local UtilityRoutines"]];
(* Read in the separate files that make up the full package. *)
(* Files in the Local directory will override files in the main directory *)
(* If the local OverlapFitting file exists, it will override the public OverlapFitting*)
If[Catch[Quiet[Check[Get["OverlapFittingLocal`"],Throw[True]],Get::noopen]],Get["OverlapFitting`"],Print["Get Local OverlapFitting"],Print["Get Local OverlapFitting"]];



End[]; (* End `Private` *)

EndPackage[];
