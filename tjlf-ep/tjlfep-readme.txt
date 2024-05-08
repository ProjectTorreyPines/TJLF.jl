This is a brief introduction to how TJLF-EP (Threads) works and how to run it.

TJLF-EP (Threads) is a version of TJLF-EP that uses the @threads macro for creating
a parallel execution of TJLF-EP. 

Threads in this version (April 16) is currently reproducing identical results to the MPI version, as it should.
It is still producing incorrect results for is_EP = 3 as MPI does as well. 

TJLF-EP (Threads) is more easily ran from the main.jl in contrast to MPI which is more easily ran from a specific directory
with the driver.jl.


Some current issues to consider:

1. The current version of Threads is taking about 10 minutes for a scan_n = 3 run. This is largely due to the fact that
the outer @threads macro loop is not running in parallel at the moment. There is an inside loop (kwscale_scan) that is also
using the macro. Running this inner loop in this way is quicker than just running the outside loop with the macro.


======================================================================================================================
|           The following information is shared between this version and the Threads version of TJLF-EP...           |                 
======================================================================================================================

TJLF-EP works in a very similar manor to the Fortran version TGLF-EP.

Here is a broad structure of the code:

                           PROCESS_IN = 5
                                  |
                                  V
main.jl/driver.jl ==> mainsub.jl ==> kwscale_scan.jl ==> ky.jl,
                                                         then tjlf_map()
                                                         then TJLF.run(), 
                                                         then TJLF.get_wavefunction()...

This txt file will be further edited later.





87001 - jl

109016 - Fortran