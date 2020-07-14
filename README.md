# kdm_kvals_code

This code creates the reference files for the correlated-k/k-distribution
method of radiative transfer for planetary atmospheres. This code covers 
the first few stesp in the radiative calculation recipe provide by
Mischna et al. (2012), doi:10.1029/2012JE004110. The code generates very 
accurate g(k) distributions for a specified gas and a range of wavenumber
for a specified range of pressures and temperatures. These g(k) 
distributions are output to files, to be read by a different code that 
merges the g(k) distribution into just a few k terms.

