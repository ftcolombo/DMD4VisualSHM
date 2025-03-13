Code implemented by F.T. Colombo for the article *Damage Identification in Structures Using Dynamic Mode Decomposition for Vibration Analysis of Low-Resolution Videos*.

# Introduction

In this work, we investigated how the Dynamic Mode Decomposition (DMD), a data-driven technique, can be interpreted in a solid mechanics problem and used to identify damage in an application of Structural Health Monitoring (SHM). The results were divided into three parts, corresponding to the numerical and experimental analysis of the motion of a beam and the implementation of a diagnostic tool.

# Numerical analysis of the motion of a clamped-free beam

The motion of a clamped-free beam was simulated according to the Euler-Bernoulli beam theory and used to create two datasets: one based on the displacement of the nodes and another based on the motion captured by a vide record.

To evaluate the results illustrated in this section, you can compile the codes in the Folder *NumericalAnalysis* in the following order:
1. analysisDMDonNodalDispl.m
2. analysisDMDonVideoMotion_undamped.m
3. analysisDMDonVideoMotion_damped.m

With DMD, we decomposed the free response of a clamped-free beam described by images of the motion into dynamic modes, which are coherent to the modes shapes seen in Structural Dynamics:

https://github.com/user-attachments/assets/7234b638-388f-42f6-bd1b-f725c32bda09

Moreover, DMD can be used to reconstruct the data in the synthetic videos used in the training. Below, we see the reconstructed motions of a beam with no damping:  

https://github.com/user-attachments/assets/a23d96bd-4068-4cdc-8560-0c0d6de4fe44

and with 4% of proportional damping in all modes:

https://github.com/user-attachments/assets/333163e1-f7dd-4050-b068-ac7579431ae1

# Experimental analysis of the motion of a clamped-free beam

Video records from the experiment done by Garrido et al.<sup>1</sup> were used to investigate the ability of DMD to work with noisy video of the motion of a flexible beam. This dataset<sup>1</sup> is needed to compile some of the implemented codes.

To evaluate the results illustrated in this section, you can compile the codes in the Folder *ExperimentalAnalysis* in the following order:
1. dmd_number_embeddings_data.m
2. dmd_number_embeddings_analysis.m
3. dmd_stabilization_data.m
4. dmd_stabilization_analysis.m
5. expvideo_dmd.m
6. expvideo_initial_guess4optdmd.m
7. expvideo_optdmd.m

Since DMD estimated conservative eigenvalues with high damping on this noisy dataset, we also investigated the algorithm *optimized DMD* proposed by Askham and Kutz<sup>2</sup>. Their MATLAB package optDMD<sup>3</sup> is required to compile some of the implemented codes.

Below, we see the motion of the beam in the videos reconstructed by DMD and optDMD (we used the eigenvalues estimated by DMD with time embeddings as the initial guess needed for optDMD):

https://github.com/user-attachments/assets/6e7fd2fc-7965-42e9-9710-640ec8a67a6c

# Diagnostic tool for SHM

To evaluate the results illustrated in this section, you can compile the codes in the Folder *DiagnosticTool* in the following order:
1. initial_modeselection4diagnostics.m
2. optdmd4diagnostics_scenarioA.m
3. diagnostics_aval_eigs_splane.m
4. diagnostics_aval_modes.m

https://github.com/user-attachments/assets/5c1815ed-0912-476a-a225-27d1856900b1

# References

1. Garrido H, Codina R, de Borbon F et al. Damage identification in beams using burst video-records during free vibration. Mechanical Systems and Signal Processing 2023; 200. URL [https://doi.org/10.1016/j.ymssp.2023.110539](https://doi.org/10.1016/j.ymssp.2023.110539)
2. Askham T and Kutz JN. Variable Projection Methods for an Optimized Dynamic Mode Decomposition. SIAM Journal on Applied Dynamical Systems 2018; 17(1): 380–416. URL [https://doi.org/10.1137/M1124176](https://doi.org/10.1137/M1124176)
3. Askham T. duqbo/optdmd: optdmd v1.0.1. Zenodo 2017. DOI 10.5281/zenodo.826433. URL  [https://github.com/duqbo/optdmd](https://github.com/duqbo/optdmd)
