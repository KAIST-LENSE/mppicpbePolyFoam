# MP-PIC-PBE_PolyFoam
This is an open-source software which is called the multiphase particle in cell coupled with the population balance equation method implemented in OpenFOAM.
The software is an extension of mppicFoam, and a tutorial is PMMA suspension polymerization which is attached together.

# Copyright Information
Copyright (C) 2021 Shin-Hyuk Kim, Jay H. Lee and Richard D. Braatz

# License
mppicPbePolyFoam is free software.
Users can redistribute it or modify it under the terms of the GNU general public license.

# Requirements
OpenFOAM 7.x is required to compile the software.

# References
Users are referred to below references for details on the MP-PIC-PBE method.

[1] S. H. Kim, R. D. Braatz, and Jay H Lee, “Multi-Phase Particle-In-Cell Coupled with Population Balance Equation (MP-PIC-PBE) Method for Multiscale Computational Fluid Dynamics Simulation,” Comput. Chem. Eng., 134, 106686, 2020.

[2] S. H. Kim, R. D. Braatz, and Jay H Lee, “Multi-scale fluid dynamics simulation based on MP-PIC-PBE method for PMMA suspension polymerization,” Comput. Chem. Eng., 152, 107391, 2021.

# Notice
When working 'wmake', some errors can occur regarding link files in the lnInclude folder at intermediate folder. In this case, replace the existing link files with the link files in openfoam5/scr/lagrangian/intermediate/lnInclude.
