# Hysteretic sticking plus CIL
Discrete element simulations of repulsive disks with a hysteretic sticking interaction with three different types of activity. This repository supports the analysis of the following manuscript:

## ABP and CIL-C
Both the ABP and CIL-C model can be compiled without any additional libraires.
```
g++ run_sd_hys_ABP.c -O3 -o MD
g++ run_sd_hys_CIL_C.c -O3 -o MD
```

## CIL-P
The CIL-P method uses a Voronoi diagram to find potential binding partners, using the Voro++ software: https://math.lbl.gov/voro++/
Compile the code using the Voro++ library.
```
g++ run_sd_hys_lang_voro_at2_ratchet_prod_poly.cpp -I path/to/voro++/voro++-0.4.6/src/ -L path/to/voro++/voro++-0.4.6/src/ -lvoro++ -O3 -o MD 
```
