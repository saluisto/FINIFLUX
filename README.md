FINIFLUX
FINIFLUX is an implicit Finite Element model that numerically solves the
steady state mass balance equation for Radon (Rn) in lotic systems.
Degasing is represented using two very popular models which originally are
from O'Connor and Dobbins (1958) and Negulescu and Rojanski (1969) with
modifications from Cartwright et al. (2011). Alternativley, FINIFLUX also can work
with degasing fluxes specified by the user. FINIFLUX is intended to
estimate groundwater fluxes into river systems as well as hyporheic
exchnage characteristics based on measured Rn concentrations. The model is
coupled to the optimasation software package PEST (Doherty 2010) for
inverse parameter estimation. FINIFLUX is intended to help scientists and
authorities that use the Rn tecnique to estimate surface-groundwater
exchange for river systems at the reach scale.
Authors: S. Frei (sven.frei@uni-bayreuth.de)
       B.S. Gilfedder (benjamin-silas.gilfedder@uni-bayreuth.de)
June 2015
June 2016: We included the posssiblity to account for surficial Rn inflow
via tributaries or drainage.
Mai 2017: We added uppwinding Finite Element correction to reduce the effect of numerical dispersion.
Hyporheic areas can now be represented using 1) exponential, 2) power law
and 3) gamma distribted RT. The entire code can be compiled using
MATLAb's application compiler. 
February 2022: fixed a bug in the degassing model 2
