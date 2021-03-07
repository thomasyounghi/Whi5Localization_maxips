# Whi5Localization_maxips
This repository contains analysis scripts for assessing G1 duration and budding interval duration in two yeast strains containing Whi5-YFP.
The two strains both have SSAcontrol cassettes, but differ in that yTY159b has a stable RFP in the cassette, while yTY160a has an RFP-degron.
The RFPdegron was expected to result in shorter G1 durations in yTY160a since the degron itself is also involved in degradation of Cln2, which is needed to progress from G1 to S phase.

## Order to run the scrips
1_CombiningData.R
2_CleaningDataforFl.R
2_IdentifyingCellsToCircle.R
3_FlProcessing.R
4_PlottingTrajectories.R
4_AssesingWhi5Thresholds.R
4_MeasuringG1duration.R
9_PlottingAgeatDox
