# TGLFEP
`TGLF` code called for EP driven local linear AE simulations

---

`TGLFEP` code is a program 
to run `TGLF` with 
energetic paritcle (EP) conveniently, 
where using the code `TGLF` (inside `GAcode`) 
with the version number 
(on May 7th, 2016) 
''TGLF stable\_r6.0.0-31-gba49'' 
on NERSC ''EDISON\_CRAY Linux x86\_64'' 
or '' CORI Linux x86\_64''.

`TGLFEP_interface.f90` contains 
the module of internal parameters
for `TGLFEP` code.

`TGLFEP_tglf_map.f90` contains 
the module of input profiles 
and the subroutine to
read plasma and geometry profiles
from input file `input.profile`
as well as the subroutine 
to map the parameters to `TGLF` 
in the required formats.

`TGLFEP_driver.f90` will read 
the control parameters 
from input file `input.TGLFEP`. 

The other subroutines contain 
different processes for 
different simulation scheme.

### Input: profiles
Four different example profiles have been tested: 
`GA-std case`, 
`DIII-D NBI case` (with GYRO inputs), 
`EPtran profiles`,
`Two EP species case`.

The existed `input.profile` is 
the example profiles for 
the DIII-D NBI case 
with 50 radii (r/a = 0.0  excluded), 
which are calculated by EPtran code. 

### Input: control parameters
The control parameters of TGLFEP code 
are read from file `input.TGLFEP`.

The existed `input.TGLFEP` has been set 
to get the EP critical beta for 
the inner 40 radii. 
The outputs are in the file `out.TGLFEP`. 
For the outside 10 radii, 
save the above output file at first 
and just change `40 SCAN_N` to `10 SCAN_N`, 
`1 IRS` to `41 IRS`, 
and `4.0 FACTOR_IN`  to `8.0 FACTOR_IN` 
in `input.TGLFEP` file. 
Then run TGLFEP again. 

`SCAN_N` is the amout of radii to calculate 
(don't have to be the same with the total number of radii, 
but should be less or equal to that). 
`IRS` is the starting radius.

`FACTOR_IN` indicates 
the upper boundary 
of the EP density scale. 
If the input EP density cannot drive AE modes, 
a larger FACTOR\_IN > 1.0 (up to 2., 10., or even 100.) 
should be used.

### EP critical beta of New cases
When you want to run a new case with TGLFEP, 
for example still the DIII-D NBI case 
but with different profiles. 
- At first, 
you should generate a new `input.profile` 
based on your plasma and geometry profiles 
in the TGLF formats. 
- Then try to run TGLFEP code and 
get the new output file. 
- Check the outputs. 
If at some radii, 
EP critical beta are not found 
(that is, `SFmin` shows NOT correctly 
in `out.TGLFEP`), 
you need to change `FACTOR_IN` in `input.TGLFEP`
to a larger value 
and try for these radii again. 

By the way, 
however it is possible that 
for some cases, 
the EP critical beta doesn't exist and
surely cannot be found 
even if with a very large `FACTOR_IN`.

### CPUs
**Recommend** to use 25x`SCAN_N` CPUs (or more) 
to run these cases. 
It usually costs about 20 mimutes 
on `Cori` (using ceil(25x`SCAN_N` /32) nodes) 
or on `Edison` (using ceil(25x`SCAN_N` /24) nodes).

---

> [More introductions on TGLF input formats](https://fusion.gat.com/theory/Tglfinput)

> The `p_prime` in `TGLF` is usually negative and should use the total pressure including the EP species.
