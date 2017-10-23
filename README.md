# Astrophysical radiative transfer

This repository hosts the source code for astrophysical radiative transfer of diffuse (wide angular
distribution) and point-source (coming from stellar sources) radiation on nested (AMR) grids. The code
was written with galaxy formation in mind but can be modified to any astrophysical problem. The radiative
transfer algorithm was described in:

* diffuse solver: Razoumov A. and Cardall C.Y., 2005, MNRAS, 362, 1413.
* point-source solver: Razoumov A. and Sommer-Larsen J., 2006, ApJ, 651, L89.

The derivatives of this code have been tested for many astrophysical problems including strong scattering
in stellar atmospheres (with good convergence).

If you need help with adapting this code to your particular astrophysical problem, please contact the
author.
