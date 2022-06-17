# Complete Observability of Nonlinear Systems
Codes for testing the complete and functional observability of nonlinear systems, based on symbolic computation of Lie derivatives.


# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".


# Usage

- `obsvmatrix` : Computes the observability matrix of a nonlinear system f(x) with a measurement function h(x) based on symbolic computations of the Lie derivatives of f(x) with respect to h(x).

- `example_rossler` : Short example for the analysis of the observability of the Rossler system, a famous chaotic system.

- `example_pendulum` : Short example for the analysis of the identifiability of the pendulum equation, an example of differential-algebraic equation. The identifiability condition can be tested based on an extended observability condition where parameters are treated as state variables. 


# References
1.  A. N. Montanari, L. A. Aguirre. Observability of Network Systems: A Critical Review of Recent Results. *Journal of Control, Automation and Electrical Systems*, **31**(6):1348–1374 (2020). DOI: 10.1007/s40313-020-00633-5.
