# TJLF
Tglf in Julia Learned from Fortran (TJLF)

Things to note: in the InputTJLF struct, WIDTH_SPECTRUM and EIGEN_SPECTRUM are new features from TGLF used to cut down eigensolver calls. KY_SPECTRUM is also an addition that saves the ky grid. FIND_WIDTH is used differently than in TGLF, if false, TJLF skips calls to tjlf_max.jl since you provide the WIDTH_SPECTRUM that would be calculated in tjlf_max.jl. FIND_EIGEN is a new flag that tells the code if you also provide a EIGEN_SPECTRUM from previous runs so that you can use eigs() instead of ggev!() in tjlf_eigensolver.jl to find the eigenvalues. This feature is used when changing gradients during flux matching a specific radial points. Changing the gradients affects the widths and eigenvalues very little allowing you to reduce calls to eigenvalue solvers which typically costs ~90% of the runtime.

Currently, Arpack.jl's eigs(), used in tjlf_LINEAR_SOLUTION.jl and tjlf_eigensolver.jl, is a fast iterative solver. In tjlf_LINEAR_SOLUTION.jl, eigs() is used to solve for the eigenvector instead of rewriting the generalized eigenvalue problem as a system of linear equations like TGLF which is potentially a little troublesome. It is also used in tjlf_eigensolver.jl if you provide a EIGEN_SPECTRUM and set FIND_WIDTH to false. eigs() does a shift and inverse iteration to solve the generalized eigenvalue problem, specifying sigma tells the solver to look for eigenvalues near sigma, it works best if you set which=:LM telling eigs() the type of eigenvalue to compute. I have tried which=:LR (LM is largest magnitude and LR is largest real part), but sometimes the solver will find "fake" eigenvalues with very large real and imaginary parts. I have found using :LM and setting sigma as the most unstable mode from the first run works.

NOTE: If you are getting:

Error: XYAUPD_Exception: Maximum number of iterations taken. All possible eigenvalues of OP has been found.
│ IPARAM(5) returns the number of wanted converged Ritz values.
│   info = 1

Make sure you are using Arpack v0.5.3 and NOT v0.5.4, the current version does not work. Arpack.jl's eigs() is NOT thread safe. I have locks in the code to keep things safe. In the future, GenericArpack.jl should provide a pure Julia version of the Arpack algorithm that is thread safe, but it is still under development and seems to be a ways off.

There are some 3D and 5D arrays where the indices are not obvious. They are specified in the function comments where they appear, but I will repeat them here:
QL_weights::Array{5} - (field, species, mode, ky, type)
    type: (particle, energy, torodial stress, parallel stress, exchange)
firstPass_eigenvalue::Array{3} - (mode, ky, type)
    type: (gamma, frequency)
QL_flux_out::Array{3} - (field, species, type)
    type: (particle, energy, torodial stress, parallel stress, exchange)
The order of the indices was to try and take advantage of Julia's column major memory usage
