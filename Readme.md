# Description for uploaded data

This repository provides additional data associated with the manuscript "Data-Driven Design of Flexible Metal-Organic Frameworks for Gas Storage (Lim et al., 2024)". </br>

There are two directories: 1) data and 2) scripts </br>

"data" directory contain the data from molecular simulations such as GCMC and DFT calculations. "scripts" directory contain the code to construct the database that we used in this work.

---

### 1. Details on each file and directory in "data" directory

#### 1.1 - 1_gpugcmc_with_zeoplusplus.csv

The .csv file that contain the coarse information of 38,562 frameworks that can be directly obtained from the name of framework (topology and components), the geometric information calculated using zeo++, and GCMC caluclations results from GPU GCMC tool [1].



#### 1.2 - 2_cds_high_mofs_with_RASPA.csv

The .csv file that contain the speciflc frameworks with "cds topology". The structures in this file satisfy certain performance threshold.</br>

Threshold: gas uptake at 65 bar > 187.2 v(STP)/v (90% of current world record (208 v(STP)/v) [2]) from both GPU GCMC and RASPA simulation tools.



#### 1.3 - 3_EV_curves

This directory contains the .txt file that record energy profiles along with volume range and associated snapshots that were obtained from the CONTCAR file from DFT calculations for each volume.</br>

There are three sub-directories and one input file for DFT calculations using VASP: 1) 0_IRMOF1, 2) 1_hypothetical, 3) 2_experimental, and 4) INCAR </br>

1. 0_IRMOF: The results of IRMOF-1 as a test case for rigid MOFs
2. 1_hypothetical: The results of 9 hypothetical MOF candidates. The visualization of components were depicted in Figure 4A. 
3. 2_experimental: The results of 9 hypothetical MOFs (MIL53_Al frameworks were conducted with two different pathway. Details are written in the manuscript.). 
4. INCAR: General INCAR files for VASP calculations.



#### 1.4 - 4_lp_np_structures

This directory contains the large pore (lp) phase and the narrow pore (np) phase frameworks for Co bdp and 9 hypothetical MOF candidates. Details on selection of volume for each phase were written in Supporting Information. Those .cif files were used for additional GCMC calculations to obtain adsorption isotherm.



#### 1.5 - References

[1]  J. Kim and B. Smit, "Efficient Monte Carlo Simulations of Gas Molecules inside Porous Materials", J. Chem. Theory Comput. 2012

[2] F. Gándara and O.M. Yaghi et al., "High Methane Storage Capacity in Aluminum Metal-Organic Frameworks", J. Am. Chem. Soc., 2014

---



---

### Citations

citation will be updated...

---

### Acknolwedgements

This work was funded by Saudi Aramco-KAIST CO2 Management Center.