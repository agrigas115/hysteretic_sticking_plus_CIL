# Hysteretic sticking plus CIL
Discrete element simulations of repulsive disks with a hysteretic sticking interaction with three different types of activity. This repository supports the analysis of the following manuscript:

## ABP and CIL-C
Both the ABP and CIL-C model can be compiled without any additional libraires.
```
g++ run_sd_hys_ABP.c -O3 -o MD
g++ run_sd_hys_CIL_C.c -O3 -o MD
```
The command line arguments for the ABP simulation are:


argv[1]: run number

argv[2]: poff

argv[3]: v0

argv[4]: phi

argv[5]: T

argv[6]: Gamma

argv[7]: seq

argv[8]: Dr


The command line arguments for the CIL_C simulation are:


argv[1]: run number

argv[2]: poff

argv[3]: v0

argv[4]: phi

argv[5]: T

argv[6]: gamma

argv[7]: seq

Units in the code are all in raw simulation units. The ABP model never generates the tensioned fluid state.
The following is an example run that generates the tensioend fluid state in the CIL-C model:
```
./MD 0 0.001 1.0 0.4 0 1000 0
```

## CIL-P
The CIL-P method uses a Voronoi diagram to find potential binding partners, using the Voro++ software: https://math.lbl.gov/voro++/
Compile the code using the Voro++ library.
```
g++ run_sd_hys_lang_voro_at2_ratchet_prod_poly.cpp -I path/to/voro++/voro++-0.4.6/src/ -L path/to/voro++/voro++-0.4.6/src/ -lvoro++ -O3 -o MD 
```
