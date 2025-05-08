(* ::Package:: *)

PacletObject[
  <|
    "Name" -> "KerrRingdown",
    "Description" -> "Kerr Ringdown Fitting Paclet",
    "Creator" -> "Leda Gao and Gregory B. Cook",
    "Version" -> "1.1.1",
    "WolframVersion" -> "12.1+",
    "PublisherID" -> None,
    "License" -> "MIT",
    "Extensions" -> {
      {
        "Kernel",
        "Root" -> "Kernel",
        "Context" -> {{"KerrRingdown`", "KerrRingdown.wl"},
        {"ReadWaveforms`", "ReadWaveforms.wl"},
        {"DataRoutines`", "DataRoutines.wl"},
        {"QNMModeData`", "QNMModeData.wl"},
        {"OverlapFitting`", "OverlapFitting.wl"},
        {"UtilityRoutines`", "UtilityRoutines.wl"}
        }
      },
      {
        "Kernel",
        "Root" -> "Local",
        "Context" -> {{"DataRoutinesLocal`", "DataRoutinesLocal.wl"},
        {"ReadWaveformsLocal`", "ReadWaveformsLocal.wl"},
        {"QNMModeDataLocal`", "QNMModeDataLocal.wl"},
        {"OverlapFittingLocal`","OverlapFittingLocal.wl"},
        {"UtilityRoutinesLocal`", "UtilityRoutinesLocal.wl"}}
      },
      {"Documentation", "Language" -> "English"},
      {"Path", "Root"->"BBH0305Compact"},
      {"Path", "Root"->"CCE0001Compact"},
      {"Path", "Root"->"QNMsdataCompact"}
    }
  |>
]
