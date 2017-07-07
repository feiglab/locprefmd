# locprefmd
local geometry protein structure refinement via MD

## 1. Installation
1.0 Prerequisites
 * CHARMM
    * http://charmm.chemistry.harvard.edu/
 * MMTSB
    * https://github.com/mmtsb/toolset/
 * MolProbity
    * https://github.com/rlabduke/MolProbity
    * Note: requires Java and php

1.1 Getting locPREFMD
 * Cloning the locPREFMD github repository
    * git clone https://github.com/feiglab/locprefmd.git

1.2 Setting environment variables
  * You have to set the following environment variables: 
    * MMTSBDIR: path to the MMTSBDIR home directory
    * CHARMMEXEC: path to the executable file of CHARMM 
    * MOLPROBITY: path to the top of the MolProbity tree
    * LOCPREFMD: path to the locPREFMD github repository after checking out

## 2. How to use locPREFMD
2.1. Prepare a input protein structure in PDB format
2.2. Run locPREFMD
  * $LOCPREFMD/scripts/locprefmd.sh [INPUT PDB] > [REFINED PDB]

## 3. Release log
  * Jul, 2017: The first release of locPREFMD

## 4. References
  * M. Feig: Local protein structure refinement via molecular dynamics simulation with locPREFMD. J. Chem. Inf. Model. (2016) 56, 1304-1312

## 5. Contact
  * feig@msu.edu

