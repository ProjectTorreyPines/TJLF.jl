# TJLF
Tglf in Julia Learned from Fortran (TJLF)

# InputTJLF

The InputTJLF struct can be populated with an InputTGLF struct, the translation is done in the tjlf_modules.jl file, but there are parameters missing between the two. The structure can also be populated with a input.tjlf file (not .gen) if you pass the file directory to the readInput() function. 

Currently, InputTJLF structs do NOT have any default values and should throw an error if it is not properly populated, this is to ensure the user is fully aware of the parameters they are using.

### New Parameters

 WIDTH_SPECTRUM and EIGEN_SPECTRUM are new features from TGLF used to cut down eigensolver calls. KY_SPECTRUM is also an addition that saves the ky grid. FIND_WIDTH is used differently than in TGLF, if false, TJLF skips calls to tjlf_max.jl since you externally provide the WIDTH_SPECTRUM that would be calculated in tjlf_max.jl. FIND_EIGEN is a new flag that tells the code if you also provide an EIGEN_SPECTRUM from previous runs so that you can use eigs() instead of ggev!() in tjlf_eigensolver.jl to find the eigenvalues. This feature is typically used when changing gradients during flux matching a specific radial points. Changing the gradients affects the widths and eigenvalues very little allowing you to reduce calls to eigenvalue solvers which typically costs ~90% of the runtime.

### Deleted Parameters

I deleted the parameters in the InputTJLF struct that are old and not used. I left comments and deleted them in commit b9a9c99b40a3ba5fcbc683eb94df4ec8b6fd5671 and commit 18f810fdf277a691b584f987781b976cccdcbb5b

# Eigenvalue Solver

Currently, Arpack.jl's eigs(), used in tjlf_LINEAR_SOLUTION.jl and tjlf_eigensolver.jl, is a fast iterative solver. In tjlf_LINEAR_SOLUTION.jl, eigs() is used to solve for the eigenvector instead of rewriting the generalized eigenvalue problem as a system of linear equations like TGLF which is potentially a little troublesome. It is also used in tjlf_eigensolver.jl if you provide a EIGEN_SPECTRUM and set FIND_WIDTH to false. eigs() does a shift and inverse iteration to solve the generalized eigenvalue problem, specifying sigma tells the solver to look for eigenvalues near sigma, it works best if you set which=:LM telling eigs() the type of eigenvalue to compute. I have tried which=:LR (LM is largest magnitude and LR is largest real part), but sometimes the solver will find "fake" eigenvalues with very large real and imaginary parts. I have found using :LM and setting sigma as the most unstable mode from the first run usually works, but more testing would be good.

```
+--------------+-----------------------------------------------+
|              |                  FirstPass()                  |
+--------------+---------+-----------------+-------------------+
|              |         | ggev!()         | eigs()            |
|              +---------+-----------------+-------------------+
|              | ggev!() |      1.141s     |                   |
|              | eigs()  |     robust,     |                   |
| SecondPass() |         | most "correct", |                   |
| (eigenvalue, |         |  thread scales? |                   |
| eigenvector) +---------+-----------------+-------------------+
|              | ggev!() |      1.091s     |                   |
|              | gsev!() |       TGLF      |                   |
|              |         |     robust,     |                   |
|              |         |  thread scales! |                   |
|              +---------+-----------------+-------------------+
|              | eigs()  |                 |     522.710 ms    |
|              |         |                 |      robust?      |
|              |         |                 |  thread scales??  |
|              |         |                 | best for 1 thread |
+--------------+---------+-----------------+-------------------+
```

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
The order of the indices was to try and take advantage of Julia's column major memory usage
