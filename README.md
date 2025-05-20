# Event-driven Monte Carlo for polydisperse hard spheres

This is an event-driven Monte Carlo (EDMC) code to perform simulations of polydisperse hard sphere in the semi-grand canonical ensemble. 

This project is adapted from a previous implementation of EDMC to deal with monodisperse particles interacting with the WCA potential (see [this release](https://github.com/FSmallenburg/EDMC) and [this reference](https://doi.org/10.1063/5.0209178) for more details).
For random number generation, it employs an implementation of Xoshiro256+ (see [this reference](https://doi.org/10.1145/3460772) and [this documentation](https://prng.di.unimi.it/)).

This code is published alongside the paper *Freezing line of polydisperse hard spheres via direct coexistence simulations*.


## Building and running

A simple Makefile is included as well as running instructions and an initial configuration for a $6%$ polydisperse system in a fluid-FCC coexistence.
All simulations parameters are defined as global parameters, some of which can be controlled with command flags (run ``make && ./md -h`` in the directory to see available flags).
User can refer to comments for more details.


## Simulation details, units, and outputs

We simulate three-dimensional systems of $N$ hard spheres in the semi-grand canonical ensemble, i.e. with fixed volume $V$, temperature $T$ and chemical potential difference $\delta\mu$.
Each particle $i$ has a position $\mathbf{r}_i$, velocity $\mathbf{v}_i$, mass $m$ (set equal for all particles), diameter $\sigma_i$ and radius $R_i = \sigma_i / 2$.
The latter is treated as a degree of freedom for the system so that particles can fluctuate in size according to interactions with the externally applied field $\delta \mu (\sigma)$.
The Hamiltonian for this system can then be written as:


$$H=\sum_i\left[\frac{\mathbf{p}_i^2}{2m} + \frac{\wp_i^2}{2M} + V(R_i)\right]$$

$$+\sum_{i\lt j}U_{ij}(\mathbf{r}_i, \mathbf{r}_j, R_i,R_j),$$

where $\mathbf{p}_i$ is the translational momentum of particle $i$, $m$ is the (translational) particle mass, $m$ is the (translational) particle mass, $\wp_i$ is the momentum associated with its radius $R_i$, $M$ is the mass associated with the radius of a particle, $V(R_i) = -\delta\mu\left(\sigma_i=2R_i\right)$ is the external field controlling the particle size, and $U_{ij}$ represents the pair interaction (which is either 0 or $+\infty$ for hard spheres).

Interactions with the continuous field $V(R_i)$ are handled with EDMC, see [this reference](https://doi.org/10.1103/PhysRevE.85.026703) for details.

<!-- Initialization of polydisperse systems with a Gaussian distribution of sizes is achieved deterministically using the inverse cummulative probability distribution function.  -->

The simulation code operates in the following units:
-  Lengths are measured in units of the mean particle size $\bar{\sigma}$.
-  Mass is measured in units of the particle mass $m$.
-  Time is measured in units of $\tau = \sqrt{\beta m \bar{\sigma}^2}$. Here, $\beta = 1/k_B T$ with $k_B$ Boltzmann's constant.
-  Energy is measured in units of $k_B T$.
-  Pressure is measured in units of $k_B T / \bar{\sigma}^3$.

<!-- The simulation code measures the pressure tensor $P$ during the simulation, and outputs it in the form of a reduced pressure $P^* = \beta P \bar{\sigma}^3$. -->
The average pressure tensor in a given time interval $[t_a, t_b]$ is measured via the virial expression

$$P_{ab} = \rho k_B T + \frac{1}{3V(t_b-t_a)} \sum  \delta p\_{i,a} r\_{ij,b}$$

where  $\rho = N/V$ is the number density with $V$ the system volume and $N$ the number of particles, and the sum in the last term is taken over all collisions in the time interval $[t_a, t_b]$.
For each collision, $\mathbf{r}_{ij}$ is the center-to-center vector connecting the two colliding particles $i$ and $j$, and $r\_{ij,b}$ is its $b$-th component ($b$ is $x$, $y$, or $z$). Similarly, $\mathbf{\delta p}_i$ is the momentum change of particle $i$ due to the collision, and $\delta p\_{i,a}$ is its $a$-th component.


## Snapshot file format

The configuration files from the simulation are written in a simple text-based format, which can contain multiple snapshots per file.
For each frame, the format consists of $N+2$ lines (with $N$ the number of particles), as follows:
- One line containing just the number of particles
- One line containing the box size, specifying the box length $L_x$, $L_y$, and $L_z$ along the three axes, separated by whitespace.
- One line per particle containing: a letter indicating particle type, three numbers indicating the real-space particle coordinates, and one number indicating the particle radius. 

The movie files that are created by the simulation code include multiple of these frames consecutively in a single text file.
Note that although the code assumes periodic boundary conditions, coordinates of particles that leave the box during the simulation will be printed as being outside of the simulation box, to allow for analysis of long-time dynamics.
Hence, any structural analysis or visualization should apply periodic boundary conditions explicitly.
This behavior is broken when the flag ``-k`` is used (equivalent to setting ``usekickback = 1``): coordinates of particles are updated so as to keep the center of mass of the crystal slab fixed during the simulation. Identification of solid-like particles is then made using ten Wolde's bond-order parameter (see [this reference](https://doi.org/10.1039/FD9960400093)).

The simulation codes can read in snapshots in this format as initial configurations (by specifying the path to the initial configuration file with the flag ``-I`` alongside the option ``-i 4`` which corresponds to ``initchoice=4``).
Periodic boundaries will be applied to the snapshot at the start of the simulation. 
Note again that the simulation code assumes that all box lengths, positions, and radii are given in units of $\bar{\sigma}$, which is the average particle diameter.
Adaptation to different configuration file formats can be done via modification of the ``loadparticles``, ``write``, and ``outputsnapshot`` functions.
