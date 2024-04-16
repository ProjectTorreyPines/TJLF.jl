This is a brief introduction to the way that TJLF-EP (MPI) works and how to run it:

TJLF-EP (MPI) is a version of TJLF-EP that utilized a message passing interface.
It is a much more direct translation of TGLF-EP and is thus similarly fast.
A single run on 3 scans with 3 processes should take about 6 minutes on a standard computer.

TJLF-EP (MPI) can be ran via a main.jl (Julia REPL) file or via the tjlfep_driver.jl file from a directory outside of the IDE.
When running this version of TJLF-EP (MPI), you must not use more than scan_n # of processes with the following command:

mpiexec -n # julia --project /PATH/TO/TJLF.jl/tjlf-ep-mpi/tjlfep_driver.jl 

where # is replaced by the number of MPI processes you want to run.
Optimally, this could then be ran on a cluster and it would be similarly fast to TGLF-EP in Fortran on the Omega cluster.



A few major issues to take note of:

This issue is still being worked out and I'm not sure why it is doing this; I have been trying to figure it out for a while. I might
be overlooking something obvious, but it seems to break down in the following places:

1. If you are running 1 more process than scan_n, (say 4 if you have scan_n = 3), the lowest of the values of ir (= 2) will seem to be "overwritten"
while the other outputted values (from out.TGLFEP) are unchanged.

2. These outputs are derived from the kwscale_scan.jl function and any discrepancies come between k = 2 and k = 3. Rounds k = 1 and k = 2 have 
the correct inputs for what would be expected.

This is why I believe there is something I'm overlooking. The kwscale_scan depends on the number of processes used for a single value of ir (np_local).
This is the ONLY thing that is different between this when running it with 1 process per ir. It should collect the values correctly with AllReduce.

For now, this current commit (April 16), will not have this issue resolved. I hope to figure out what is actually wrong with the code in the next few weeks.


How to run TJLF-EP in the MPI format:

As previously mentioned, there are two ways to do this:

    - main.jl
    - tjlfep_driver.jl

The first of these is mostly for testing purposes. You can run a specific file in a chosen directory there from the TJLF.jl directory itself. When running from here, you will find
your output files loose inside of the TJLF.jl directory rather than a remote location. main.jl is ran with REPL.

The latter of these is better for containing results. You can run a specific input.TGLFEP and input.MTGLF (see below) from a directory on your computer. This is where the command
"mpiexec -n # julia --project /PATH/TO/TJLF.jl/tjlf-ep-mpi/tjlfep_driver.jl" is used. 

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

