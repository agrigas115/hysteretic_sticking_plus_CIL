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

argv[2]: poff - probability of unbinding

argv[3]: v0 - swim velocity

argv[4]: phi - packing fraction

argv[5]: T - temperature in Langevin thermostat

argv[6]: gamma - friction in Langevin thermostat

argv[7]: seq - int for radii set in seq/seq_%s.txt

argv[8]: Dr - activity vector rotational diffusion coefficient


The command line arguments for the CIL_C simulation are:


argv[1]: run number

argv[2]: poff - probability of unbinding

argv[3]: v0 - swim velocity

argv[4]: phi - packing fraction

argv[5]: T - temperature in Langevin thermostat

argv[6]: gamma - friction in Langevin thermostat

argv[7]: seq - int for radii set in seq/seq_%s.txt

Units in the code are all in raw simulation units. The ABP model never generates the tensioned fluid state.
The following is an example run that generates the tensioend fluid state in the CIL-C model:
```
./MD 0 0.001 1.0 0.4 0 1000 0
```
Output is saved in two files in directories traj_data and bonded_data. The trajectory files contain the x and y positions of each particle, as well as the x and y components of each particles activity vector, respectively arranged in four columns. Each frame is separated by a line that lists the box size four times. The bonded data contains the particle indices are bonded pairs, separated by tabs. Each new line corresponds to the next frame. Example frames are provided. Every 100,000 steps, checkpoint files are written to chk to allow automatic restarting. Deleted previous checkpoint files if starting a new run. The seed is set by run number and sequence.

## CIL-P
The CIL-P method uses a Voronoi diagram to find potential binding partners, using the Voro++ software: https://math.lbl.gov/voro++/
Compile the code using the Voro++ library.
```
g++ run_sd_hys_lang_voro_at2_ratchet_prod_poly.cpp -I path/to/voro++/voro++-0.4.6/src/ -L path/to/voro++/voro++-0.4.6/src/ -lvoro++ -O3 -o MD 
```
The command line arguments for the CIL_P simulation are:

argv[1]: seq - int for radii set in seq/seq_%s.txt

argv[2]: run number

argv[3]: pon - probability of making a new adhesion

argv[4]: poff - probability of cutting an exisiting adhesion

argv[5]: phi - packing fraction

argv[6]: T - temperature in Langevin thermostat

argv[7]: gamma - friction in Langevin thermostat

argv[8]: theta - cutoff in adhesion angle to determine if a new adhesion is possible

argv[9]: rr - ratchet rate of pulling in adhesion rest length
