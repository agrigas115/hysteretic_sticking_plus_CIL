# Hysteretic sticking plus CIL
Discrete element simulations of repulsive disks with a hysteretic sticking interaction and contact inhibition of locomotion (CIL). This repository supports two manuscripts:

1. *Sparse mesenchymal cell networks as a fluid under tension*: https://www.biorxiv.org/content/10.64898/2025.12.07.692626v1. Code at the top level.
2. *Fluids under tension: dynamics, mechanics, and structure of sparse active adhesive networks*: https://arxiv.org/abs/2611.31838. Code in `fluids_under_tension/`.

If you find this code useful, please cite the manuscripts. Please contact me if you have any questions.

## Sparse mesenchymal cell networks as a fluid under tension
All simulations can be visualized with their corresponding plot.py function.

### ABP and CIL-C
Both the ABP and CIL-C model can be compiled without any additional libraries.
```
g++ run_sd_hys_lang_ABP.c -O3 -o MD
g++ run_sd_hys_lang_CIL_C.c -O3 -o MD
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
The following is an example run that generates the tensioned fluid state in the CIL-C model:
```
./MD 0 0.001 1.0 0.4 0 1000 0
```
Output is saved in two files in directories traj_data and bonded_data. The trajectory files contain the x and y positions of each particle, as well as the x and y components of each particle's activity vector, respectively arranged in four columns. Each frame is separated by a line that lists the box size four times. The bonded data contains the particle indices of bonded pairs, separated by tabs. Each new line corresponds to the next frame. Example frames are provided. Every 100,000 steps, checkpoint files are written to chk to allow automatic restarting. Delete previous checkpoint files if starting a new run. The seed is set by run number and sequence.

### CIL-P
The CIL-P method uses a Voronoi diagram to find potential binding partners, using the Voro++ software: https://math.lbl.gov/voro++/
Compile the code using the Voro++ library.
```
g++ run_sd_hys_lang_CIL_P.cpp -I path/to/voro++/voro++-0.4.6/src/ -L path/to/voro++/voro++-0.4.6/src/ -lvoro++ -O3 -o MD
```
The command line arguments for the CIL_P simulation are:

argv[1]: seq - int for radii set in seq/seq_%s.txt

argv[2]: run number

argv[3]: pon - probability of making a new adhesion

argv[4]: poff - probability of cutting an existing adhesion

argv[5]: phi - packing fraction

argv[6]: T - temperature in Langevin thermostat

argv[7]: gamma - friction in Langevin thermostat

argv[8]: theta - cutoff in adhesion angle to determine if a new adhesion is possible

argv[9]: rr - ratchet rate of pulling in adhesion rest length

The following is an example run that generates the tensioned fluid state in the CIL-P model:
```
./MD 0 0 0.01 0.001 0.4 0 1000 2.0 0.1
```

## Fluids under tension: dynamics, mechanics, and structure of sparse active adhesive networks
`fluids_under_tension/` holds:
- `cil_bd.c`: both hysteretic models, crawling (`--model C`) and pulling (`--model P`), in one code. It is the same model as above, rewritten for speed, still single core: neighbour and cell lists, bond lifetimes drawn at birth, binding as a continuous-time process, and its own radical Voronoi construction (no Voro++).
- `attractive_cil.c`: the attractive particle model with CIL-C, which falls into an absorbing state.
- `trajectory.py`: reader for the `traj.bin` snapshots of `cil_bd`.

Neither code needs libraries:
```
cd fluids_under_tension
make
```

### Units
Overdamped dynamics at T = 0 with drag 1 and spring constant 1, integrated with forward Euler. Rates are per unit time: koff (bond lifetime tau_b = 1/koff), and for CIL-P kon and kr. The codes above (dt = 0.1, gamma = 1000) convert as t = t_old/20, koff = -20 ln(1 - poff)/100, kon = -20 ln(1 - pon)/100, kr = 20 rr/100. The paper uses f0 = 1 (CIL-C), kon = 0.021, kr = 0.02, theta_m = 2 (CIL-P), and phi = 0.4. These are the defaults.

### Bulk runs
Random positions at phi, relaxed by FIRE, then t_equil of dynamics and t_run of measurement. The paper uses t_equil = 5 tau_b and t_run = 50 tau_b. For example, at koff = 0.001:
```
./cil_bd --model C --seed 0 --N 1024 --koff 0.001 --t_equil 5000 --t_run 50000
./cil_bd --model P --seed 0 --N 1024 --koff 0.001 --t_equil 5000 --t_run 50000
```
The seed sets the radii (uniform in [0.5, 0.7)), the positions and all random events.

### Body-force probes (moduli)
A probe starts from the equilibrated state a bulk run saves (`equilibrated.txt`). Each cell feels a force F0 sin[q(y_i - Y)] along x (transverse, `--vF0`) or V0 q sin[q(y_i - Y)] along y (longitudinal, `--V0`), with q = 2 pi/L and Y the centre of mass. The paper uses a peak force of 0.05 q (F0 = 0.05 q, V0 = 0.05), runs each probe at +F0 and -F0 from the same state, and measures for 20 tau_b. For example, N = 256 (q = 0.232, so F0 = 0.01162) at koff = 0.001:
```
D=data/phi_0.4_f0_1_koff_0.001_N_256_dt_0.02/seed_0
./cil_bd --model C --koff 0.001 --start $D/equilibrated.txt --probe 1 --vF0 0.01162 --t_strain 20000 --v_split 1 --v_log 20 --v_every 250 --traj 0
./cil_bd --model C --koff 0.001 --start $D/equilibrated.txt --probe 1 --vF0 -0.01162 --t_strain 20000 --v_split 1 --v_log 20 --v_every 250 --traj 0
./cil_bd --model C --koff 0.001 --start $D/equilibrated.txt --probe 1 --V0 0.05 --t_strain 20000 --v_split 1 --v_log 20 --v_every 250 --traj 0
./cil_bd --model C --koff 0.001 --start $D/equilibrated.txt --probe 1 --V0 -0.05 --t_strain 20000 --v_split 1 --v_log 20 --v_every 250 --traj 0
```
Each writes `vf.txt` in `$D/probe_<T|L>_n1_<F|V><amp>/`, the mode displacement and force modes against time, from which G(t), K(t) and eta follow.

### Output
Runs write to `data/phi_<phi>_f0_<f0>_koff_<koff>_N_<N>_dt_<dt>/seed_<seed>/` (CIL-C) or `data/phi_<phi>_kon_<kon>_kr_<kr>_koff_<koff>_N_<N>_dt_<dt>/seed_<seed>/` (CIL-P), with each stage (`bulk`, a probe, or a strain step) in its own directory:
- `isf.txt`, `chi4.txt`, `msd.txt`: self-intermediate scattering function at q = 2 pi sqrt(rho), four-point susceptibility, and mean-squared displacement (centre of mass removed), on a log-spaced lag grid.
- `response.txt`: pair stress (Pxx, Pyy, Pxy) and mean coordination against time.
- `traj.bin`: snapshots every 0.1 tau_b and log-spaced bursts, with bonds (and, for CIL-P, the rest lengths of bonds still pulling in). Read with `trajectory.py`.
- `vf.txt`: probe mode velocity, displacement and density modes.
- `log.txt`, `params.txt`, `radii.txt`, `equilibrated.txt`, `endpoint.txt`.

Runs checkpoint every 30 minutes and on SIGTERM, and resume exactly from `checkpoint.bin`; a finished run exits at once. Delete the run directory to start over. `./cil_bd --help` lists every option, including step strains (`--shear`, `--dilation`) and extending a finished run (`--extend 1`).

To read the snapshots in Python:
```
from trajectory import read_frames, bonds_to_pairs
for frame in read_frames('data/.../bulk/traj.bin'):
    pairs = bonds_to_pairs(frame['bonds'], len(frame['radii']))
```
Each frame gives unwrapped positions, radii, the box size, time and the bonds; `bonds_to_pairs` turns the bonds into cell pairs (i, j). CIL-P frames also give `rest_lengths`.

### Attractive particle model
Each cell has a repulsive core of radius R_i and an attractive shell of radius lambda R_i, so a pair is repulsive while the cores overlap and attractive within the shells, with no hysteresis. Each cell crawls with force v0 away from the sum of its contact directions (CIL-C). The shell packing fraction phi_a sets lambda = sqrt(phi_a/phi). The paper uses phi_a = 0.9, a crawling force v0 = 2.5 (lambda - 1) <sigma_ij> (v0 = 1.5 at phi = 0.4), 20 seeds, t_run = 1e6, and varies the core packing fraction phi and the time step dt. For example:
```
./attractive_cil --seed 1 --phi 0.4 --phi_a 0.9 --v0 1.5 --dt 0.005 --t_run 1e6 --log_every 20000 --traj_every 100000
```
Output goes to `data/phi_<phi>_phia_<phi_a>_v0_<v0>_N_<N>_dt_<dt>/seed_<seed>/`: `log.txt` (contacts, speeds, bond turnover and pair stress against time), `window.txt` (fraction of active cells per window), and `traj.txt` (frames of x, y and the crawling direction of each cell, each frame starting with a line holding the box size). `./attractive_cil --help` lists every option.
