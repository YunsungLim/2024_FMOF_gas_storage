# Description for uploaded data

This repository provides additional data associated with the manuscript **"Data-Driven Design of Flexible Metal-Organic Frameworks for Gas Storage (Lim et al., 2024)"**. </br>

There are two directories: 1) data and 2) scripts </br>

"data" directory contain the data from molecular simulations such as GCMC and DFT calculations. "scripts" directory contain the code to construct the database that we used in this work.

---



## 1. Details on "data" directory :floppy_disk:

### 1.1 - 1_gpugcmc_with_zeoplusplus.csv

- The .csv file that contain the coarse information of 38,562 frameworks that can be directly obtained from the name of framework (topology and components), the geometric information calculated using zeo++, and GCMC caluclations results from GPU GCMC tool [1].



### 1.2 - 2_cds_high_mofs_with_RASPA.csv

- The .csv file that contain the speciflc frameworks with "cds topology". The structures in this file satisfy certain performance threshold.

- Threshold: gas uptake at 65 bar > 187.2 v(STP)/v (90% of current world record (208 v(STP)/v) [2]) from both GPU GCMC and RASPA simulation tools.



### 1.3 - 3_EV_curves

- This directory contains the .txt file that record energy profiles along with volume range and associated snapshots that were obtained from the CONTCAR file from DFT calculations for each volume.

- There are three sub-directories and one input file for DFT calculations using VASP: 1) 0_IRMOF1, 2) 1_hypothetical, 3) 2_experimental, and 4) INCAR 

1. 0_IRMOF: The results of IRMOF-1 as a test case for rigid MOFs
2. 1_hypothetical: The results of 9 hypothetical MOF candidates. The visualization of components were depicted in Figure 4A. 
3. 2_experimental: The results of 9 hypothetical MOFs (MIL53_Al frameworks were conducted with two different pathway. Details are written in the manuscript.). 
4. INCAR: General INCAR files for VASP calculations.



### 1.4 - 4_lp_np_structures

- This directory contains the large pore (lp) phase and the narrow pore (np) phase frameworks for Co bdp and 9 hypothetical MOF candidates. Details on selection of volume for each phase were written in Supporting Information. Those .cif files were used for additional GCMC calculations to obtain adsorption isotherm.



### 1.5 - References

[1]  J. Kim and B. Smit, "Efficient Monte Carlo Simulations of Gas Molecules inside Porous Materials", J. Chem. Theory Comput. 2012

[2] F. Gándara and O.M. Yaghi et al., "High Methane Storage Capacity in Aluminum Metal-Organic Frameworks", J. Am. Chem. Soc., 2014

---



## 2. Details on 'scripts' directory :computer:

Every codes were tested in python 3.8 and pormake 0.2.0.</br>

Highly recommend not to run following scripts in the original environment (in user interface) with pormake because there is some modification to pormake!!



### 2.1 - 1_modified_codes_for_pormake

- 3 files (builder.py, framework.py, and locator.py) of original pormake library were revised in this work. Copy this 3 files into the installed pormake directory.

  - builder.py: build() function in the Builder() class was revised

  - framework.py: cutoff for neighbor_list() was changed as 12.0 from 6.0 (line 139)

  - Locator.py: locate_with_permutation() function of Locator() class was revised




### 2.2 - 2_notebook

#### - "rodmof" directory

- Additional files for rodmof generation.

- The topology files (.cgd file format) and rod building block files (.xyz file format) are provied in the assets directory. The name of asymmetric edge building blocks that we exclude in this work is also provided in the assets directory as .csv file format.

#### - Identify_topo.ipynb

- Jupyter notebook to extract 'cds', 'dia', and 'mog' topology that we used in this manuscript.

- The raw .cif files for extraction were provided in the 'assets' directory (not 'assets' directory in the 'rodmof directory'). As written in the manuscript, we calculated the coordination sequence from the raw cif file using the duplicating metal atom of rod metal node building blocks as pseudo-bonding site. 

- Details on scratching the topology is provided in the manuscript and comments in the notebook file.

#### - Single_generation_example.ipynb

- Jupyter notebook to provide example for generation of single framework.

- The representative wine-rack MOFs, Co bdp, MIL53(Al), and MIL118, were generated as test case.

#### - Bulk_generation.ipynb

- Jupyter notebook to generation of frameworks in bulk scale (Database construction).

- Total 39,005 frameworks can be constructed, and 7 mog frameworks were eliminated because geometry optimization failure using Materials Studio. The number and the exact shape of frameworks can be slightly different from every different environment (version of python libraries, geometry optimization process, or etc).

- This workflow or notebook file can also be used as general bulk generation using PORMAKE with slight modification.



### 2.3 - 3_assets

- .cif file of 3 representative framework (Co bdp, MIL53(Al), and MIL118) for running jupyter notebook named Identify_topo.ipynb in notebook directory.
- Perl scripts file for geometry optimzation of frameworks using Materials Studio

---



## 3. Citations :page_with_curl:

[Y. Lim, B. Kim, and J. Kim*, Data-Driven Design of Flexible Metal-Organic Frameworks for Gas Storage, *Chem. Mater.*, 2024](https://pubs.acs.org/doi/full/10.1021/acs.chemmater.4c00398)

---



## 4. Acknolwedgements :muscle:

This work was funded by Saudi Aramco-KAIST CO2 Management Center.