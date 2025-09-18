## Markovian entropy balance

This repository simulates Langevin trajectories under periodic feedback control and computes several thermodynamical quantities: work, heat, entropy, housekeeping entropy, mutual information, transfer entropy and unavaiable information.
Several parameters of the protocol (f, $\Delta$<sub>tm</sub>, $\Delta$<sub>x</sub>, V<sub>0</sub>) can be customized to explore different regimes. Comparations of different bounds for work based on information quantities are plotted.

# Distribution of files


* <ins>Fortran folder</ins>: This folder contains files to simulate the Langevin equation with feedback protocol that can includes measurements errors. *Feedback_trajectory.f90* is the main file. *protocol_Wm.f90*, *probability.f90*, *binary_fun.f90*, and *ran3.f* files are subroutines.
* <ins>Matlab folder</ins>: This folder contains files to compute the thermodynamical quantities and compare the different bounds based on information quantities. 
