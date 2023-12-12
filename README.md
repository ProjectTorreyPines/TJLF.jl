# TJLF
Tglf in Julia Learned from Fortran (TJLF)

# InputTJLF

The InputTJLF struct can be populated with an InputTGLF struct, the translation is done in the tjlf_modules.jl file, but there are parameters missing between the two. The structure can also be populated with a input.tglf file (not .gen) if you pass the file directory to the readInput() function. 

Currently, InputTJLF structs do NOT have any default values (maybe like one or two exceptions) and should throw an error if it is not properly populated, this is to ensure the user is fully aware of the parameters they are using.

## New Parameters

InputTJLF mostly follows the input.tglf file variables, but there are a couple new things

FIND_WIDTH::Bool - this is in the input.tglf file, but is used different than in the Fortran, I have repurposed this variable to be a flag to tell the code whether or not we are providing an external width spectrum to avoid extra calls to tjlf_max.jl that solves for width (see [Eigenvalue Solver](#eigenvalue-solver))

WIDTH_SPECTRUM::Vector - external width spectrum used if FIND_WIDTH is False. If FIND_WIDTH is True, this Vector will be populated with widths for each ky grid point. The idea is that you can run with FIND_WIDTH set as True, it will automatically save the width spectrum, you can then change the gradients of that particular struct, and then change the FIND_WIDTH flag to False. Or, you can easily save this value to other InputTJLF struct, not sure why you would necessarily, but the information is here.

FIND_EIGEN::Bool - this flag tells the code which solver to use when solving for the eigenvalues, if unsure, FIND_EIGEN as True is more robust (see [Eigenvalue Solver](#eigenvalue-solver)). Should ONLY be False if FIND_EIGEN is also False because it would be suspicious if you are confident on the eigenvalues enough to use and not confident on the width values (it would also probably break tjlf_max.jl), there should be an assertion to check this.

EIGEN_SPECTRUM::Vector - external eigen spectrum used on the first pass if FIND_EIGEN is False. This value is set to the sigma parameter of the eigs() solver and acts as an initial guess for the iterative solver. Currently, this is a vector that automatically stores the most unstable eigenvalues from the first pass, but it might be smart to save the eigenvalues for each mode and rework how eigs works in LS.jl

EIGEN_SPECTRUM2::Vector - similar to above, but save the most unstable eigenvalues from the second pass, after the spectral shift since this shifts the eigenvalues

SMALL::Float - value used to rewrite the eigenvalue problem as a system of linear equations like the Fortran does (default 1e-13)

## Deleted Parameters

I deleted the parameters in the InputTJLF struct that are old and not used. I left comments and deleted them in commit [b9a9c99](https://github.com/ProjectTorreyPines/TJLF.jl/commit/b9a9c99b40a3ba5fcbc683eb94df4ec8b6fd5671) and commit [18f810f](https://github.com/ProjectTorreyPines/TJLF.jl/commit/18f810fdf277a691b584f987781b976cccdcbb5b).

Deleted parameters: USE_TRANSPORT_MODEL,GEOMETRY_FLAG,B_MODEL_SA,FT_MODEL_SA,VPAR_SHEAR_MODEL,WRITE_WAVEFUNCTION_FLAG,VTS_SHEAR,VNS_SHEAR,VEXB,RMIN_SA,RMAJ_SA,Q_SA,SHAT_SA,ALPHA_SA,XWELL_SA,THETA0_SA,NN_MAX_ERROR

# Eigenvalue Solver

Currently, Arpack.jl's eigs(), can be used in tjlf_LINEAR_SOLUTION.jl and tjlf_eigensolver.jl as a fast iterative solver. In tjlf_LS.jl, eigs() is used to solve for the eigenvector instead of rewriting the generalized eigenvalue problem as a system of linear equations like TGLF, which is potentially a little troublesome. It is also used in tjlf_eigensolver.jl if you provide a EIGEN_SPECTRUM and set FIND_EIGEN to False. eigs() does a shift and inverse iteration to solve the generalized eigenvalue problem, specifying sigma gives the shift, by setting which=:LM you tell eigs() to compute eigenvalues around sigma. I have tried which=:LR (LM is largest magnitude and LR is largest real part), but sometimes the solver will find "fake" eigenvalues with very large real and imaginary parts. I have found using :LM and setting sigma as the most unstable mode has worked, but more testing would be good.

Table of the three different combinations of the solvers used in TJLF
```
+--------------+-----------------------------------------------+
|              |                  FirstPass()                  |
+--------------+---------+-----------------+-------------------+
|              |         | ggev!()         | eigs()            |
|              +---------+-----------------+-------------------+
|              | ggev!() |      1.141s     |                   |
|              | eigs()  |     robust,     |                   |
|              |         | most "correct", |                   |
|              |         |  thread scales? |                   |
|              +---------+-----------------+-------------------+
| SecondPass() | ggev!() |      1.091s     |                   |
| (eigenvalue, | gsev!() |       TGLF      |                   |
| eigenvector) |         |     robust,     |                   |
|              |         |  thread scales! |                   |
|              +---------+-----------------+-------------------+
|              | eigs()  |                 |     522.710 ms    |
|              |         |                 |      robust?      |
|              |         |                 |  thread scales??  |
|              |         |                 | best for 1 thread |
+--------------+---------+-----------------+-------------------+
```
To run top left, set the InputTJLF SMALL parameter = 0.0, set FIND_EIGEN = False<br>
To run mid left, use the default InputTJLF SMALL parameter = 1e-13<br>
To run bot right, set the InputTJLF SMALL parameter = 0.0, set FIND_EIGEN = True<br>

# Arpack.jl

NOTE, If you are getting:<br>
<pre>
Error: XYAUPD_Exception: Maximum number of iterations taken. All possible eigenvalues of OP has been found.<br>
│ IPARAM(5) returns the number of wanted converged Ritz values.<br>
│   info = 1</pre>

Make sure you are using Arpack v0.5.3 and NOT v0.5.4, the current version does not work. Also, Arpack.jl's eigs() is NOT thread safe. I have locks in the code to keep things safe. In the future, GenericArpack.jl should provide a pure Julia version of the Arpack algorithm that is thread safe, but it is still under development and seems to be a ways off.

# Indices of Arrays

There are some 3D and 5D arrays where the indices are not obvious. They are specified in the function comments where they appear, but I will repeat them here:<br>
<pre>QL_weights::Array{5} - [field, species, mode, ky, type]<br>
    type: (particle, energy, torodial stress, parallel stress, exchange)</pre><br>
<pre>firstPass_eigenvalue::Array{3} - [mode, ky, type]<br>
    type: (gamma, frequency)</pre><br>
<pre>QL_flux_out::Array{3} - [field, species, type]<br>
    type: (particle, energy, torodial stress, parallel stress, exchange)</pre><br>
The order of the indices try to take advantage of Julia's column major memory usage
