# Virginia Tech Common Research Wind Tunnel (VT CRWT) inflow velocity generator

This repository contains codes corresponding to the VT CRWT project.

The purpose of the software is to generate a two-dimensional velocity field by fitting experimental data. The profile can be generated and adjusted by running and editing the profile_gen_2D.m MATLAB script.


### TODO

- [] improving curve fitting near the wall
- [] taking asymemtry into account

### References
Further details are provided in the work of Althaf et al. “Numerical modeling and tunnel specific considerations for CFD model development of low-speed wind tunnels”, AIAA AVIATION Forum, San Diego CA, USA, 2023​
https://doi.org/10.2514/6.2023-3980​

**Please consider citing the aforementioned study when you use this repository for your work.**

### Dependencies
The code relies on two open-source repositories which are redistributed based on the terms of the corresponding licenses:
1; Suraj Shankar (2023). 2D Poisson equation (https://www.mathworks.com/matlabcentral/fileexchange/38090-2d-poisson-equation), MATLAB Central File Exchange. Retrieved October 9, 2023.
2; Audrey Cheong (2023). Euclidean distance between two point clouds (https://www.mathworks.com/matlabcentral/fileexchange/59377-euclidean-distance-between-two-point-clouds), MATLAB Central File Exchange. Retrieved October 9, 2023.

Repository 2 was repurposed as a function to obtain the poisson_solver_2D.m

### Contributors
Affiliations during development are presented
Carola Rovira Sala - Cranfield University
Tamás István Józsa - Cranfield University (tamas[DOT]jozsa[AT]cranfield[DOT]ac[DOT]uk)
Máté (Matt) Szőke - Virginia Tech

CRS developed the software under TIJ's supervision. Experimental data (TravLog_v1_Summary50_MAT.mat) and a preliminary profile esimate (parabola_fit.mat) were provided by MSz.
