In this work, we investigated how the Dynamic Mode Decomposition (DMD), a data-driven technique, can be interpreted in a solid mechanics problem and used as a diagnostic tool for Structural Health Monitoring (SHM).  The results were divided into three parts, corresponding to the numerical and experimental analysis of the motion of a beam and the implementation of a diagnostic tool.

## Numerical analysis of the motion of a clamped-free beam

The motion of a clamped-free beam was simulated according to the Euler-Bernoulli beam theory and used to create two datasets: one based on the displacement of the nodes and another based on the motion captured by a vide record.

To evaluate the results illustrated in this section, you can compile the codes in the Folder *NumericalAnalysis* in the following order:
1. analysisDMDonNodalDispl.m
2. analysisDMDonVideoMotion_undamped.m
3. analysisDMDonVideoMotion_damped.m

With this numerical analysis, we can decompose the motion described in the synthetic videos into dynamic modes, which are coherent to the modes shapes seen in Structural Dynamics for this type of beam:

https://github.com/user-attachments/assets/6c087826-458d-4fa8-871b-8d27ef7585d5

Moreover, DMD can be used to reconstruct the data in the synthetic videos used in the training. Below, we see the reconstructed motions of a beam without damping and with 4% of proportional damping in all modes: 

https://github.com/user-attachments/assets/771881da-3883-4fd5-9c13-aa8b53d6b2ac

https://github.com/user-attachments/assets/b1dd6070-517f-4b39-8be1-594def039161

## Experimental analysis of the motion of a clamped-free beam

Video records from the experiment done by x were used to investigate the ability of DMD to work with noisy video of the motion of a flexible beam. This dataset is needed to compile some of the implemented codes.

To evaluate the results illustrated in this section, you can compile the codes in the Folder *ExperimentalAnalysis* in the following order:
1. dmd_number_embeddings_data.m
2. dmd_number_embeddings_analysis.m
3. dmd_stabilization_data.m
4. dmd_stabilization_analysis.m
5. expvideo_dmd.m
6. expvideo_initial_guess4optdmd.m
7. expvideo_optdmd.m

Since DMD estimated conservative eigenvalues with high damping on this noisy dataset, we also investigated the algorithm *optimized DMD* proposed by Askham and Kutz. Their MATLAB package is required to compile some of the implemented codes.

Below, we see the beam motion in the videos reconstructed by DMD and optDMD, where we used the eigenvalues estimated by DMD as the initial guess needed for optDMD:

## Diagnostic tool for SHM

To evaluate the results illustrated in this section, you can compile the codes in the Folder *DiagnosticTool* in the following order:
1. initial_modeselection4diagnostics.m
2. optdmd4diagnostics_scenarioA.m
3. diagnostics_aval_eigs_splane.m
4. diagnostics_aval_modes.m

