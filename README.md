# Dimensionless Numbers Edupack


## Instructions For Use

**If you do not have MATLAB installed**, download *installer_DNEdupack.exe* to install the application and run it.

**If you already have MATLAB installed**, download the whole repository in a *.zip* format and extract it. All the files should be in the same directory. Run *DNEdupack.mlapp* in MATLAB.


## Design Philosophy

This prototype presents a software toolkit that facilitates the learning and appreciation of dimensionless numbers in fluid mechanics through the use of interactive simulations. In our project, wherein we are designing an aerosol extraction system, we find dimensionless numbers such as the Stokes number and the Ohnesorge number indispensable in quantifying distinct flow regimes and predicting aerosol behaviour. Students are first introduced to dimensional analysis in MECH0005: Introduction to Thermodynamics and Fluid Mechanics. However, as they are presented with numerous dimensionless numbers accompanied only by brief explanations of their physical meaning, students often fail to develop a robust intuition or appreciation of their significance.

This edupack remedies that deficiency by providing interactive simulations and supplementary theoretical reading material [1], thereby enabling students to cultivate both conceptual understanding and practical insight into the governing dimensionless parameters of fluid mechanics.

By opening the edupack, students will be presented with a menu where they can select the dimensionless number which most interests them. 

<img width="508" height="380"  alt="image" src="https://github.com/user-attachments/assets/4c7c02eb-4949-45c6-8945-e2c2ae1cb1e8" /> 

Upon selection, they will be given a brief explanation as to the physical meaning of the dimensionless number and what flow physics it encompasses.

<img width="508" height="380"   alt="image" src="https://github.com/user-attachments/assets/fd26b4b9-b7eb-4ae7-8493-b4eb42186965" />

Students can then proceed to the simulation to observe the effect of the dimensionless number of fluid flow and vary parameters to observe the response.

<img width="508" height="380"   alt="image" src="https://github.com/user-attachments/assets/90cb400f-53ff-4e5f-86dd-33d629482342" /> 

The simulations resolve the governing equations behind fluid mechanics with an emphasis on computational speed, thus necessitating the use of reduced order models where appropriate. The simulations were numerically resolved with a combination of public available [2,3] and proprietary code using the Finite Difference Method, Finite Volume Method or the Spectral Method where most appropriate. For the purposes of this edupack, the targeted run speed was on the order of millisecond to seconds to permit interactive or near-real-time visualisation of flow phenomena on standard consumer hardware. Whereas concessions to numerical fidelity have been made, the results of the simulations are still qualitatively accurate in showing the variation across the pertinent dimensionless numbers if not exactly quantitatively accurate. 


## References

[1] Lau NCK, Klettner CA. Vorticity generation, transport, and annihilation under steep solitary waves. Physics of Fluids. 2025; 37(9):097173. https://doi.org/10.1063/5.0287404 

[2] MathWorks-Teaching-Resources. Computational-Fluid-Dynamics [Internet]. 2025 [cited 2025 Oct 13]. Available from: https://github.com/MathWorks-Teaching-Resources/Computational-Fluid-Dynamics/tree/main

[3] Ricciardi A.. AndersonCFD_Chapter10_Solution (version 1.0) [Internet]. 2021 [cited 2025 Oct 13]. Available from: https://uk.mathworks.com/matlabcentral/fileexchange/95203-andersoncfd_chapter10_solution



