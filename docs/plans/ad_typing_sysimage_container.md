# TJLF: type-stability first, then sysimage, then container

## Context

TJLF's ForwardDiff path costs ~20 minutes the first time it is used
(`TJLF/test/runtests.jl:20` says so in as many words). Having just shipped TJLFEP
containers and sysimages, the question is whether to clone that pipeline for TJLF.

**Answer: worthwhile, but staged — and neither a container nor a sysimage is the fix.**
The lead item is type stability, per your read that proper variable typing should help the
AD solve. Exploration says that is right, and localizes it.

### What the 20 minutes is, and what each option does to it

The 20 min is *compile* work, cached per
`(package version × dependency manifest × Julia version × CPU target)`, and re-paid
whenever any of those change. Evidence it keeps being re-paid: this depot holds **9 distinct
pkgimage slugs** each for `TJLF` and `TJLFForwardDiffExt` under `v1.11` at ~85–91 MB apiece
(~1.6 GB), plus separate `v1.10` and `v1.12` trees. NERSC's default Julia module has since
moved to **1.12.1**, which invalidates all of it again.

| | Effect on the 20 min |
|---|---|
| **Type stability** | Reduces real work — in the conversion/actor layer at compile time, and in *every* AD evaluation at runtime. The only item that also makes the solve itself faster. |
| **Sysimage** | Does not reduce the compile; it *is* the compile, snapshotted. Pays off across many workers/jobs in one campaign. |
| **Container** | Ships that snapshot to people who would otherwise rebuild it. The only way to guarantee an outside user never pays it, because it freezes exactly one environment. |

### Where the typing debt actually is

`InputTJLF{T}` is **already clean** — the `Union{...,Missing}` purge happened, and
`tjlf_modules.jl:236-241` documents the recipe and the reason ("a `Union{Missing,...}` here
poisons inference across all of TJLF and forces heap-boxing of millions of scalars"). The
hot structs (`Ave`, `AveH`, `AveWG`, `OutputGeometry`, …) are cleanly parametric
`Array{T,3}` / `Vector{T}`. Only 5 `ismissing` calls remain in all of `src/`.

The remaining debt is concentrated in four places, all on the Dual path:

1. **`InputTGLF{T}` — 142 of 179 fields are `Union{T,Missing}`**
   (`tjlf_modules.jl:3-234`). FUSE constructs this **with `T = Dual`**:
   `tglf_actor.jl:93` builds `InputTGLF{D}`, then `:104`/`:106` reads it field-by-field.
   Note FUSE imports the name from TurbulentTransport, but
   `TurbulentTransport.jl:7` re-exports TJLF's — there is one definition, and it lives in
   TJLF. This is the same fix `InputTJLF` already received.

2. **`minimal_scalar_copy` (`tjlf_modules.jl:513-521`) is fully dynamic and runs inside the
   ky loop.** It loops `for f in fieldnames(typeof(inputs))` doing `setfield!`/`getfield`
   with a **runtime** Symbol — every access infers to `Any`. It is called at
   `tjlf_TRANSPORT_MODEL.jl:42,53,63,73,83`, i.e. inside `Threads.@threads for ky_index`,
   across two passes — roughly `179 fields × 2·nky` boxed field operations per solve, and
   *per AD evaluation*.

3. **`update_input_tjlf!` (`tjlf_modules.jl:473-510`)** has the same dynamic-Symbol pattern,
   including `getfield(input_tglf, Symbol("ZS_", i))` which allocates a Symbol per species
   per field.

4. **FUSE's actor field is Union-typed**:
   `input_tglfs::Union{Vector{InputTGLF{D}},Vector{InputTJLF{D}}}` (`tglf_actor.jl:31`), so
   the `TJLF.run_tjlf(actor.input_tglfs)` call at `:189` dispatches through a Union.

### The performance target that makes this measurable

The reference configuration is the flux-matcher's `:simple_trust` algorithm, which is the
one gradient-based path: `flux_matcher_actor.jl:221` resolves `jacobian_method` to
`:forward_ad` for `:simple_trust` (and `:finite_diff` for everything else), and `:322` maps
it to `NonlinearSolve.SimpleTrustRegion(; autodiff)`. Reference cases are the DIII-D L-mode
and H-mode cells of `FuseExamples/fluxmatcher.ipynb` with `act.ActorTGLF.model = :TJLF`.

**The acceptance criterion: a Dual solve should cost about 2× a Float64 solve** (versus
`N+1` solves for a finite-difference Jacobian). TJLF's AD is IFT-based
(`TJLFForwardDiffExt.jl:32,89,200,258`): the eigensystem is computed **once** in Float64
LAPACK, independent of chunk size, and the ~2× floor is concrete — the non-Hermitian
`_herm_eigen` path does a second `eigen(Af')` for the left vectors (`:108`), and
`_standard_eigenvalues_via_solve` adds two `lu` factorizations (`:287,300`).

Be precise about what *does* scale with `N`: the derivative propagation is six explicit
`for k in 1:np` loops (`:49,149,224,275,290,304`), each a Float64 BLAS product or
triangular solve per partial. That is cheap relative to an eigensolve, so the ratio should
sit near 2 at small `N` and grow **slowly and linearly** — that growth is correct behaviour,
not a defect. What the benchmark is looking for is a ratio well above that line, or one that
grows super-linearly, and the plausible causes are enumerable:

- boxing from the `Union{T,Missing}` fields and dynamic-Symbol field ops (1c/1d below),
- a fall-through to generic `eigen` on Dual inputs (1f), which the extension's own comment
  at `:245-257` calls "very slow",
- Dual arithmetic dominating in some region the IFT does not cover.

That ratio turns "proper typing should help the AD solve" into a number that either holds
or points at which of the three is responsible.

### The chunk size is currently unpinned

`flux_matcher_actor.jl:252-253` constructs `ADTypes.AutoForwardDiff()` with **no
`chunksize`**, so ForwardDiff derives it from the number of unknowns
(`channels × length(rho_transport)`). That number moves with the case: the default
`rho_transport = 0.25:0.1:0.85` (`:15`) gives 7 radii, while the notebook's benchmark cell
uses `0.1:0.05:0.85` with `evolve_rotation`/`evolve_densities = :flux_match` — 16 radii and
4 channels. Unknowns swing from ~14 to ~64, and `ForwardDiff.pickchunksize` maps those onto
different chunk sizes (~7 vs ~11).

That has two costs: per-evaluation work varies unpredictably, and **each distinct chunk size
is another full specialization of the 14.4k-LOC solve** to compile and to bake. Pinning
`AutoForwardDiff(; chunksize=…)` is a one-line change that collapses the specialization
count and makes the bake list finite and knowable. Phase 1a measures the runtime trade
before choosing the value.

### The honest split

Typing and specialization-baking address **different halves**, and it is worth being clear
before spending effort:

- Typing fixes the conversion/actor layer and **every AD evaluation's** boxing and
  allocation. This is the "help the AD solve" win, it is real, and the 2× ratio measures it.
- It does **not** remove the bulk of the first-call cost, which is instantiating the
  ~14.4k-LOC eigensolver stack for each distinct `Dual{Tag,Float64,N}`. Those structs are
  already concretely parametric, so that cost is inherent *per chunk size* — which is what
  pinning the chunk size and baking a known list address instead.

Which is why Phase 1 starts by **measuring the apportionment** rather than assuming it.

### Two coverage gaps that would make an image bake the wrong thing

Independently of typing, building a sysimage or container today would snapshot an
incomplete workload:

- **Only chunk size `N=1` is baked.** `TJLFForwardDiffExt.jl:341-348` traces
  `ForwardDiff.derivative`. Real consumers use other `N`, and each re-specializes the whole
  solve: `TJLFEP.gamma_input_sensitivities` uses `N = 3*NS+1` → 10/13/16 for NS=3/4/5
  (`tjlfep_ad_extensions.jl:1930-1942`); `gamma_grad` uses `N = length(vars)` (`:311-320`);
  TJLF's own `runtests_ad.jl:113` uses `N=2`; FUSE's chunk is set by the flux-matcher.
- **Only the scalar entry point is baked.** FUSE calls the **Vector** method —
  `tglf_actor.jl:189`, `run_tjlf(actor.input_tglfs)` → `FluxSolution{D}`.
  `run_tjlf(::Vector{InputTJLF{Dual}})` (`run_tjlf.jl:138`, `Threads.@threads`) is a
  separate specialization that is never traced.

And one thing to diagnose rather than assume: two of the nine `TJLFForwardDiffExt` caches
are **72 KB** instead of ~91 MB (one `TJLF` cache is 2 MB instead of 85 MB) — the workload
emitted almost no code there. `PrecompileTools` swallows workload errors, so a
silently-failing workload is indistinguishable from a fast one.

**User decisions:** optimize for FUSE `use_ad`, external collaborators/CI, and NERSC batch
campaigns; deliver all three phases, staged; bake only the handful of chunk sizes actually
in use.

---

## Phase 1 — Type stability + workload coverage

No node hours; login node plus at most one small CPU job. This is the phase that carries
the value.

### 1a. Benchmark and apportion (before changing anything)

**Runtime — the 2× test.** Build the harness on the pattern already sitting (commented out)
in `FuseExamples/fluxmatcher.ipynb` cell 2, which does exactly this shape of `@btime`
breakdown for the `:QLNN` path. Retarget it at `:TJLF` with `algorithm = :simple_trust`,
on the DIII-D L-mode and H-mode cases:

```julia
using BenchmarkTools
inp64 = actor.actor_ct.actor_turb.input_tglfs          # Vector{InputTJLF{Float64}}
t64  = @belapsed TJLF.run_tjlf($inp64)                  # one finite-diff evaluation
tdual = @belapsed TJLF.run_tjlf($inpD)                  # same solve, InputTJLF{Dual{…,N}}
@show tdual / t64                                       # target ≈ 2
```

Sweep `N` over the chunk sizes the two cases actually produce. Report `tdual/t64` **and**
`@ballocated` for both — allocation growth is the fingerprint of boxing, and it will move
before wall time does. Do the same for the scalar and `Vector` entry points separately,
since only the scalar one is currently precompiled.

Then localize whatever gap the ratio shows:

- `JET.@report_opt` on `run_tjlf(::InputTJLF{Dual})` and on `InputTJLF{D}(::InputTGLF{D})`;
  `@code_warntype` on `minimal_scalar_copy`.
- `Profile.@profile` (or `ProfileView`) on a Dual solve to confirm where the excess sits —
  in particular whether any Dual eigen call reaches the generic `eigen` fallback.
- Record a full `:simple_trust` flux-match wall time for L-mode and H-mode as the
  end-to-end before/after number.

**Compile — the apportionment.** In a throwaway depot so nothing warm masks the cost:

```bash
export JULIA_DEPOT_PATH=$PSCRATCH/adprobe/.julia   # fresh, disposable
```

- `SnoopCompileCore.@snoopi_deep` on a first `ForwardDiff.gradient`, then
  `flamegraph`/`accumulate_by_source`: how much of the 20 min is the eigensolver stack
  instantiating per `N`, versus the Union-splitting conversion layer?
- Time `Pkg.precompile()` for TJLF alone and with `ForwardDiff` present; then time
  first-call after a warm cache. This is the precompile-vs-JIT split that sets the value of
  Phases 2 and 3.
- Record the chunk sizes both reference cases select, as the input to the bake list and to
  the pinning decision in 1e.

### 1b. Diagnose the undersized caches

Reproduce a 72 KB `TJLFForwardDiffExt` cache and surface the swallowed error (temporarily
wrap the workload body in `try … catch e; @error e; end`, or run it outside
`@compile_workload`). Candidates: `pkgdir(TJLF)` returning `nothing`,
`precompile/sample_input.tglf` absent from an installed copy, or `Pkg.test`'s
`--check-bounds=yes` disabling pkgimages. If the workload silently fails in a common
environment, that alone explains a recurring 20 minutes.

### 1c. Make `InputTGLF{T}` concretely typed

`TJLF/src/tjlf_modules.jl:3-234`. Apply exactly the recipe already documented at `:236-241`
for `InputTJLF`: concrete field types, sentinel defaults (`T(NaN)` for floats, `0` for ints,
`false` for bools, `T[]` for vectors), with "unset" still guarded by the NaN sentinel rather
than `Missing`.

This is the largest single edit in the plan (142 fields) but it is mechanical and the
precedent is in the same file. Sequence it carefully:

- `readInput` / `tjlf_read_input.jl` must populate sentinels rather than leave `missing`.
- The `ismissing(v) ||` guard in `update_input_tjlf!:480` becomes a NaN check (or drops out).
- Only 5 `ismissing` sites exist in `src/`; audit each. One is already dead code, not a
  conversion candidate: the SHAPE_* loop at `update_input_tjlf!:501-505` calls `ismissing`
  on `InputTJLF` fields that are now concrete with `zero(T)` defaults, so it is always
  false — delete it (the intended 0.0 default is already the constructor default).
- Check downstream users of `InputTGLF` that may rely on `missing`: `TurbulentTransport`
  (re-exporter), FUSE's `tglf_actor.jl` (`:112-114` copies fields from
  `par.custom_input_files`), and TJLFEP.
- The AD conversion helpers already skip missing values defensively
  (`TJLFForwardDiffExt.jl:327`, `runtests_ad.jl`, `TJLFEP._to_dual_input`); those branches
  become dead and should go.

### 1d. Make the field copies static

Replace the dynamic-Symbol loops in `minimal_scalar_copy` (`:513-521`) and
`update_input_tjlf!` (`:473-510`) with `@generated` functions that unroll to compile-time
`QuoteNode` field symbols, so every access is typed and no Dual is boxed:

```julia
@generated function minimal_scalar_copy(inputs::InputTJLF{T}) where {T<:Real}
    assigns = [:(setfield!(out, $(QuoteNode(f)), getfield(inputs, $(QuoteNode(f)))))
               for f in fieldnames(InputTJLF)]
    quote
        out = InputTJLF{T}(inputs.NS, length(inputs.KY_SPECTRUM))
        $(assigns...)
        return out
    end
end
```

This is the change most likely to speed up the AD **solve** itself, since it sits inside the
threaded ky loop and runs on every evaluation. Two things to verify while doing it: the
current version aliases the parent's `Vector` fields into the copy (`setfield!` shares the
reference) — preserve that behaviour deliberately rather than by accident, and confirm the
`@generated` form does not measurably inflate per-`N` compile time (it unrolls ~179 typed
stores; expected to be well below what it saves).

### 1e. Extend the AD workload to the shapes actually used

`TJLF/ext/TJLFForwardDiffExt.jl:341`, `TJLF/src/TJLF.jl:117`.

- Run each target chunk size through **both** the scalar and the `Vector` entry point.
- Reuse the existing `_convert_input_tjlf` (`:323`) rather than adding a fourth copy of that
  helper — it is currently duplicated in `runtests_ad.jl:10`, the README, and TJLFEP.
- Gate the chunk list with **Preferences.jl**, not an environment variable: preferences
  participate in the cache slug, so changing the list correctly invalidates the cache,
  whereas an env var would silently produce a mismatched image. Adds `Preferences` to
  `[deps]`; if that is unwanted, hardcode the list and drop the override.
- **Default `[1, 2]`** — cheap, covers `derivative` plus small gradients. Larger sizes are
  opt-in, since each additional `N` costs ~90 MB of cache and real minutes.
- The Phase-2/3 build environments set the preference to the campaign list, e.g.
  `[1, 2, 10, 13]` plus the flux-matcher chunk from 1a.

**Pin the flux-matcher chunk size.** Change `ADTypes.AutoForwardDiff()` at
`flux_matcher_actor.jl:253` to pass an explicit `chunksize`, chosen from the 1a sweep as the
best runtime trade across the L-mode and H-mode unknown counts. This is a one-line edit that
makes the set of `Dual{Tag,Float64,N}` specializations finite and known instead of a
function of `rho_transport` — which is what lets the bake list in Phase 2/3 actually cover
production. Keep it a `Switch`/`Entry` so a user can override it, and confirm against the
1a numbers that pinning does not cost meaningful wall time on either case.

Keep `precompile/sample_input.tglf` (3-species, EM, SAT2, `NBASIS_MAX=6`). `SAT_RULE` is a
runtime value, not a type parameter — `runtests_ad.jl:132-133` already records that rules
1/2/3 compile quickly once the IFT eigen dispatch exists — so there is no need to bake
multiple SAT rules.

### 1f. Cheap adjacent check

Confirm every Dual eigen call site reaches a specialized wrapper and none falls through to
generic `eigen`; the extension's own comment at `:245-257` notes the fallback is "very slow
on Dual inputs". Sites: `tjlf_matrix.jl:430,442,473,485` and
`tjlf_eigensolver.jl:2737-2757`.

### 1g. Housekeeping

`Project.toml` has **no `julia` compat entry** while `Manifest.toml` records 1.11.7 and
NERSC now defaults to 1.12.1. Add an entry for the versions actually tested; do not claim
1.12 without running the suite there.

No downstream compat bumps are needed: FUSE (`TJLF = "1.2.4"`), TurbulentTransport
(`"1.2.3"`) and TJLFEP (`"1.2.12, 1.3"`) are caret ranges that already admit 1.3.x.

### Gate 1 → 2

If first-AD-use drops to seconds and precompile to a few minutes, the sysimage's marginal
value falls sharply. Reassess before spending node hours.

---

## Phase 2 — TJLF sysimage

New `TJLF/build/sysimage/`, mirroring `TJLFEP/build/sysimage/` (clean, committed, measured):

| New file | Adapted from |
|---|---|
| `setup_registry_env.jl` | same name — swap the package set to `TJLF`, `ForwardDiff`, `PackageCompiler` (+ `CUDA` pinned `5.8.5` for GPU) |
| `build_cpu_sysimage.jl` | same name → `create_sysimage([:TJLF, :ForwardDiff]; …)` |
| `build_gpu_sysimage.jl` | `build_gpu_sysimage_fileonly.jl` → `create_sysimage([:CUDA, :TJLF, :ForwardDiff]; …)` |
| `precompile_cpu_workload.jl` | same name — Float64 scalar + Vector `run_tjlf`, then AD at **each** target chunk size through **both** entry points |
| `batch_build_{cpu,gpu}_sysimage.sh` | the TJLFEP batch scripts |

Carry the TJLFEP lessons across verbatim:

- **`cpu_target = ENV["JULIA_CPU_TARGET"]`, not `PackageCompiler.default_app_cpu_target()`.**
  The in-repo TJLFEP builders still use the latter, which traces the build host only
  (Milan); the container builder
  (`deploy/container/gpu-sysimage/build_sysimage_container.jl`) fixes this and says so.
  Start TJLF on the fixed version.
- Keep workload bodies inside a `let`, never top-level `const` — a baked `Main.GACODE`
  collides with task scripts that define their own.
- Batch script order: `cudatoolkit/12.9` **then** `julia/1.11.7`;
  `JULIA_CUDA_USE_COMPAT=false`; `umask 007`.
- A sysimage does **not** embed JLL artifacts and re-reads package *sources* by absolute
  path. Anything shared must be baked against the **CFS depot**
  (`TJLFEP_DEPOT=$CFS/depot`), never `$PSCRATCH` or `$HOME`. Never `rsync -a` into the share
  (it preserves the wrong group).
- Set the Phase-1 chunk-size preference in the generated build env.

Accounts default to TJLFEP's (`-A m3739` CPU, `-A m3739_g` GPU, `-q regular` — m3739 lost
premium QOS); one-line change if TJLF should be charged elsewhere.

**Expected cost**, from TJLFEP's measured builds (CPU ~14 min/1.2 GB, GPU ~28 min):
~10–20 min and ~0.25 GPU-node-hours per image; size ~600 MB–1 GB (base Julia + ~85 MB
Float64 TJLF + ~90 MB per baked chunk size).

### Gate 2 → 3

Build the container only if there is a real standalone or external consumer.

---

## Phase 3 — TJLF container

Clone `TJLFEP/deploy/container/` into `TJLF/deploy/container/`: `Containerfile`,
`install_tjlf_container.jl`, `build.sh`, `publish_ghcr.sh`, `acceptance.sh`,
`install_tjlf_container_nersc.sh`, `test_slurm_{cpu,gpu}.sbatch`, `gpu-sysimage/`, and
`.github/workflows/container.yml`.

**Caveat to settle first.** TJLFEP is an *application* (file-based `input.TGLFEP` runs), so
a container is its natural unit. TJLF is a *library*: a collaborator who wants it as a
dependency of their own project cannot consume an image. To make this a real deliverable
rather than a warm depot in a box, ship a `/usr/local/bin/tjlf` CLI taking an `input.tglf`
and writing fluxes (plus a `--grad` mode exercising the AD path), mirroring TJLFEP's
launcher including its read-only-depot fallback and `--sysimage` probe. Without that CLI,
Phase 2 plus a documented CFS sysimage already covers the NERSC and FUSE audiences, and
Phase 3 buys only CI reproducibility.

Deltas from the TJLFEP pipeline:

- **One CI variant**, not two — TJLFEP's lean/imas split exists because `TJLFEPIMASExt`
  needs IMAS/GACODE/TurbulentTransport; TJLF's weakdeps are just CUDA and ForwardDiff and
  both belong in the image. Tags: `:<v>`, `:latest`, `:<v>-gpu`.
- Publish to `ghcr.io/projecttorreypines/tjlf`. **The first push creates the package
  private** — an org owner must flip it public before unauthenticated pulls work.
- Keep the CUDA 12.9 `LocalPreferences.toml` pin written **before** any CUDA package is
  added, then `install_all_artifacts(include_lazy=true)`, so the image works offline on
  compute nodes. `local = "false"` (artifact) for the container vs `local = "true"` (module
  toolkit) for bare metal.
- Keep the FuseRegistry SSH→HTTPS rewrite, the `ARG TJLF_VERSION` cache-bust, the REPL
  warm-up layer, and the CI free-disk step (~10 GB; the CUDA artifact is ~3 GB).
- Resolve TJLF with an explicit `Pkg.PackageSpec(name="TJLF", version=…)` — TJLFEP learned
  the expensive way that an unversioned `Pkg.add` can silently resolve back to a version
  that does not even precompile.
- `:latest` must never point at a `-gpu` tag (`publish_ghcr.sh` enforces this).
- `-gpu` is baked manually on an A100 from the **pulled CI image**, not a local build — the
  `.so` embeds depot paths and package versions and must match the published base.

---

## Verification

**Phase 1** — against a fresh disposable depot:
- **The headline gate: `@belapsed` ratio of a Dual solve to a Float64 solve is ~2 at small
  `N` and grows only slowly and linearly with `N`** (the O(N) IFT propagation loops), at the
  flux-matcher's chunk size on both the DIII-D L-mode and H-mode cases, with `@ballocated`
  for the Dual solve dropping sharply. If the ratio sits well above that line or grows
  super-linearly after 1c/1d, the remaining cause is in 1f (generic `eigen` fallback) or in
  a non-IFT region — report which, rather than declaring the phase done.
- End-to-end: a full `:simple_trust` flux-match on both cases, wall time vs the 1a baseline.
- `JET.@report_opt` on the Dual path is clean where it was not; `@code_warntype` on
  `minimal_scalar_copy` shows concrete types throughout.
- Compile: `Pkg.precompile` time before/after, and first `ForwardDiff.gradient` at **each**
  baked `N` through **both** entry points, completing in seconds. `--trace-compile` AD
  entries should largely disappear.
- Every `TJLFForwardDiffExt` cache is ~90 MB-class, never 72 KB.
- **Correctness is the gate on the `InputTGLF` retyping**:
  `cd TJLF && julia --project=test test/runtests.jl` in full — the 14 `tglf_regression`
  cases must still match Fortran `out.tglf.gbflux` to `atol=1e-2`, and `test_core`/`test_EM`
  golden outputs to `RTOL 5e-4`.
- Downstream smoke: a FUSE `ActorTGLF` run on both the `InputTGLF` and `InputTJLF` branches
  (`tglf_actor.jl:65/67`), plus one TJLFEP `run_gacode_scan_task` with `solver=:ad`, since
  both touch the retyped struct.
- Re-run one FUSE flux-matching case with `use_ad=true`; first iteration should no longer
  stall.

**Phase 2** — `sbatch build/sysimage/batch_build_cpu_sysimage.sh`, then under the image:
- `Base.JLOptions().image_file` ends with the expected `.so`; full regression suite passes.
- AD at each baked `N` is fast on first call (wall-time bound, as TJLFEP's GPU acceptance
  does).
- Portability: runs on a non-Milan node — exactly what `default_app_cpu_target()` gets wrong.
- Shareability: from a fresh empty primary depot stacked in front of the CFS depot, and
  `find $CFS/... ! -group m3739` comes back empty.

**Phase 3**:
- `acceptance.sh` on a login node: version identity, sysimage loaded, writable-depot
  fallback, offline (`--network none`) solve, `CUDA.functional()==false` degrading to CPU,
  registries world-readable, plus a TJLF-specific check that gradients at the baked `N`
  return without a compile stall.
- CI smoke on `ubuntu-24.04` (300 min timeout): version assert + offline regression case.
- `sbatch deploy/container/test_slurm_gpu.sbatch` on an A100 for `-gpu`.
- Anonymous `podman-hpc pull` / `apptainer pull` once the package is public.

---

## Phase 1 results (2026-08-26, login node, 8 threads, `sample_input.tglf`)

Benchmark: `benchmark/ad_ratio.jl` (best of 5, scalar entry / Vector entry of 4).

| case | scalar best | ratio | scalar alloc | vector best | ratio |
|---|---|---|---|---|---|
| Float64 | 0.81 s | 1.00 | 311 MB | 3.66 s | 1.00 |
| Dual N=1 | 3.16 s | 3.88 | 1.87 GB | 11.5 s | 3.13 |
| Dual N=2 | 3.92 s | 4.83 | 2.75 GB | 13.3 s | 3.63 |
| Dual N=7 | 7.99 s | 9.82 | 7.18 GB | 25.9 s | 7.09 |
| Dual N=11 | 11.1 s | 13.67 | 10.7 GB | 33.1 s | 9.03 |

- **The ratio sits at ~3.9 at N=1 (above the ~2 IFT floor) and grows ~linearly at
  ~0.9×/partial — before AND after 1c/1d (baseline N=1 was 3.93).** 1f confirmed no
  generic-`eigen` fallback, so per this plan's own triage the excess is Dual arithmetic
  in non-IFT regions (matrix assembly / geometry / FLR carry partials at O(N) per scalar
  op; allocation grows ~0.9 GB/partial as the fingerprint). The InputTJLF solve kernel
  was already type-stable; 1c/1d pay off in the conversion/actor layer and compile
  coverage, not in the solve kernel. At N≈11 an AD evaluation ≈ N+1 finite-diff
  evaluations in wall time; AD wins on accuracy/robustness, and on time only at small N.
  Follow-up lever if the ratio matters: push the IFT boundary outward (freeze partials
  through matrix assembly, differentiate the assembled matrices).
- **Tags matter.** First call at an unbaked *tag* costs ~40–50 s per chunk size even for
  baked `N` (inference is reused; native code is tag-specific). The workload bake
  mitigates but cannot fully cover a consumer's first call unless the consumer pins
  `AutoForwardDiff(chunksize=…, tag=…)` to match — the FUSE pin reduces the surface to
  one warm-up per session per chunk size.
- **1b resolved: the 72 KB caches are flag variants, not swallowed errors.** The 69 MB
  `.ji`-without-`.so` entries are `--check-bounds=yes` (Pkg.test) builds — that flag
  disables native caching entirely; the 71 KB `.so` entries match a `--compile=min`
  tooling process (editor language server). The workload itself runs fine.
- Precompile cost with the default `[1, 2]` chunk list × both entry points: TJLF 81 s,
  TJLFForwardDiffExt 603 s, ext cache 324 MB (was 87 MB with the old N=1-scalar-only
  workload).
- Correctness gate: full `test/runtests.jl` passes (14 tglf_regression vs Fortran, SAT,
  core/EM golden, eigen, batched SI, ForwardDiff AD 23/23 in 9m42s).

## Recommended stopping points

- **Phase 1 alone** is worthwhile regardless of the rest — it is the only phase that reduces
  total work rather than relocating it, and the only one that speeds up the AD solve itself.
  Within it: 1a is cheap and decides everything downstream; 1d (static field copies) is the
  cheapest change with the most direct effect on per-evaluation AD cost; the chunk-size pin
  in 1e is one line and makes the whole bake story tractable; 1c (`InputTGLF` retyping) is
  the largest edit and needs the regression suite as its gate.
- **Phase 2** pays for itself on multi-worker campaigns and for m3739 collaborators via CFS.
  Note a sysimage is invalidated by the same environment churn that produced nine cache
  slugs — pinning one environment per campaign matters as much as the image.
- **Phase 3** is worth it if the `tjlf` CLI makes the image a genuine standalone tool, or if
  reproducible external/CI distribution is a goal in itself. As a pure latency fix it is the
  weakest of the three.
