module TJLFForwardDiffExt

using TJLF
using ForwardDiff
using LinearAlgebra
import LinearAlgebra.LAPACK.ggev!

# ─────────────────────────────────────────────────────────────────────────────
# AD-compatible eigen for ForwardDiff.Dual element types.
#
# LAPACK only supports Float32/64/ComplexF32/ComplexF64, so for Dual-number
# elements we compute the Float64 eigensystem with LAPACK and propagate the Dual
# perturbation analytically via the implicit function theorem (IFT):
#
#   ∂λᵢ/∂p  =  lᵢᴴ (∂A/∂p) rᵢ                         (simple eigenvalue)
#   ∂rᵢ/∂p  =  Σ_{j≠i}  [lⱼᴴ (∂A/∂p) rᵢ / (λᵢ - λⱼ)] rⱼ
#
# IMPORTANT (no type piracy): these methods are defined on TJLF-OWNED functions
# (`TJLF._sym_eigen`, `TJLF._herm_eigen`, `TJLF._generalized_eigenvalues`), NOT
# on `LinearAlgebra.eigen`. TJLF's eigen call sites route through those owned
# wrappers, so AD specialization here stays local to TJLF and never overrides
# `eigen` for other packages.
# ─────────────────────────────────────────────────────────────────────────────

const _DEGEN_THRESHOLD = 1e-12

# ── Real symmetric matrix with Dual entries ──────────────────────────────────
# Matches:  modwd!  →  _sym_eigen(Symmetric(ave.wdh)) / _sym_eigen(Symmetric(ave.wdg))
# Returns a NamedTuple supporting both `wh, ah = _sym_eigen(wdh)` (positional)
# and `eig.values` / `eig.vectors` (field access).
function TJLF._sym_eigen(A::LinearAlgebra.Symmetric{D,M}) where {D <: ForwardDiff.Dual, M <: AbstractMatrix{D}}
    np  = ForwardDiff.npartials(D)
    Tag = ForwardDiff.tagtype(D)

    # Float64 value matrix (Symmetric wrapper preserved)
    Af = LinearAlgebra.Symmetric(map(ForwardDiff.value, parent(A)))
    ef = eigen(Af)                    # LAPACK – works for Float64 Symmetric
    λf = ef.values                    # Vector{Float64}
    V  = ef.vectors                   # Matrix{Float64}
    n  = length(λf)
    pA = parent(A)

    # Precompute B_k = Vᵀ (∂A/∂pₖ) V for each partial direction k
    # B_k[i,i] = ∂λᵢ/∂pₖ   and   B_k[j,i]/(λᵢ−λⱼ) = S_k[j,i]  (j≠i)
    dλ = zeros(n, np)
    dV = zeros(n, n, np)

    for k in 1:np
        dAk = LinearAlgebra.Symmetric(map(x -> ForwardDiff.partials(x, k), pA))
        Bk = V' * (dAk * V)

        # eigenvalue derivatives (diagonal of Bk)
        for i in 1:n
            dλ[i, k] = Bk[i, i]
        end

        # eigenvector derivatives via mixing-matrix S
        Sk = zeros(n, n)
        for i in 1:n, j in 1:n
            if i != j
                gap = λf[i] - λf[j]
                if abs(gap) > _DEGEN_THRESHOLD
                    Sk[j, i] = Bk[j, i] / gap
                end
            end
        end
        dV[:, :, k] .= V * Sk    # dV_k = V * S_k
    end

    # Construct Dual eigenvalues
    λ = map(1:n) do i
        ForwardDiff.Dual{Tag}(λf[i], ntuple(k -> dλ[i, k], Val(np))...)
    end

    # Construct Dual eigenvectors
    Vd = Matrix{D}(undef, n, n)
    for j in 1:n, i in 1:n
        Vd[i, j] = ForwardDiff.Dual{Tag}(V[i, j], ntuple(k -> dV[i, j, k], Val(np))...)
    end

    return (values = λ, vectors = Vd)
end

# ── Complex matrix with Dual entries ─────────────────────────────────────────
# Matches:  modkpar!  →  _herm_eigen(im .* ave.kpar) / _herm_eigen(im .* ave.kpar_eff[is,:,:])
#           eigensolver  →  _herm_eigen(B \ A)
# eigenvalues may be complex (e.g. pure-imaginary for anti-Hermitian matrices)
function TJLF._herm_eigen(A::AbstractMatrix{Complex{D}}) where {D <: ForwardDiff.Dual}
    np  = ForwardDiff.npartials(D)
    Tag = ForwardDiff.tagtype(D)

    # Float64 value matrix
    Af = map(a -> Complex{Float64}(ForwardDiff.value(real(a)), ForwardDiff.value(imag(a))), A)

    # Use Hermitian-aware solver when applicable (better-conditioned eigenvectors)
    herm = ishermitian(Af)
    ef = herm ? eigen(Hermitian(Af)) : eigen(Af)
    λf = ef.values    # Vector{ComplexF64} or Vector{Float64}
    V  = ef.vectors   # Matrix{ComplexF64}
    n  = length(λf)

    # For non-Hermitian matrices we need LEFT eigenvectors for correct
    # derivatives.  For Hermitian / normal matrices, left = right.
    if herm
        L = V   # Hermitian ⟹ orthonormal ⟹ left = right
    else
        ef_adj = eigen(Af')
        # eigen(A') gives eigenvalues conj(λ) and right-eigenvectors of A'
        # which are the LEFT eigenvectors of A.
        L_raw = ef_adj.vectors
        λ_adj = ef_adj.values          # ≈ conj.(λf)

        # Match left eigenvectors to right eigenvectors by eigenvalue
        perm = zeros(Int, n)
        used = falses(n)
        for i in 1:n
            best_j = 0
            best_d = Inf
            for j in 1:n
                if !used[j]
                    d = abs(conj(λ_adj[j]) - λf[i])
                    if d < best_d
                        best_d = d
                        best_j = j
                    end
                end
            end
            perm[i] = best_j
            used[best_j] = true
        end
        L = L_raw[:, perm]

        # Bi-orthogonal normalisation: scale so that lᵢᴴ rᵢ = 1
        for i in 1:n
            s = dot(L[:, i], V[:, i])
            if abs(s) > 1e-30
                L[:, i] ./= conj(s)
            end
        end
    end

    # Precompute derivatives for all partial directions
    dλ_re = zeros(n, np)
    dλ_im = zeros(n, np)
    dV_re = zeros(n, n, np)
    dV_im = zeros(n, n, np)

    for k in 1:np
        dAk = map(a -> Complex{Float64}(ForwardDiff.partials(real(a), k),
                                        ForwardDiff.partials(imag(a), k)), A)
        # Bk = Lᴴ dA R  (reduces to Vᴴ dA V for Hermitian)
        Bk = L' * (dAk * V)

        # eigenvalue derivatives: ∂λᵢ = Bk[i,i]
        for i in 1:n
            dλ_re[i, k] = real(Bk[i, i])
            dλ_im[i, k] = imag(Bk[i, i])
        end

        # eigenvector derivatives via mixing-matrix S
        # For Hermitian:   dV = V * S,  S[j,i] = B[j,i]/(λᵢ−λⱼ)
        # For non-Hermitian: dR = R * S,  S[j,i] = B[j,i]/(λᵢ−λⱼ)
        Sk = zeros(ComplexF64, n, n)
        for i in 1:n, j in 1:n
            if i != j
                gap = λf[i] - λf[j]
                if abs(gap) > _DEGEN_THRESHOLD
                    Sk[j, i] = Bk[j, i] / gap
                end
            end
        end
        dVk = V * Sk
        dV_re[:, :, k] .= real.(dVk)
        dV_im[:, :, k] .= imag.(dVk)
    end

    # Construct Dual eigenvalues
    λ = map(1:n) do i
        Complex(ForwardDiff.Dual{Tag}(real(λf[i]), ntuple(k -> dλ_re[i, k], Val(np))...),
                ForwardDiff.Dual{Tag}(imag(λf[i]), ntuple(k -> dλ_im[i, k], Val(np))...))
    end

    # Construct Dual eigenvectors
    Vd = Matrix{Complex{D}}(undef, n, n)
    for j in 1:n, i in 1:n
        Vd[i, j] = Complex(ForwardDiff.Dual{Tag}(real(V[i, j]), ntuple(k -> dV_re[i, j, k], Val(np))...),
                           ForwardDiff.Dual{Tag}(imag(V[i, j]), ntuple(k -> dV_im[i, j, k], Val(np))...))
    end

    return (values = λ, vectors = Vd)
end

# ── Generalized complex eigenproblem A r = λ B r with Dual entries ───────────
# Matches:  eigensolver  →  _generalized_eigenvalues(amat, bmat)  (non-transport
# model / EP single-ky path). Returns ONLY the eigenvalues (Vector{Complex{Dual}})
# to match the owned `_generalized_eigenvalues` contract.
#
#   ∂λᵢ/∂p = lᵢᴴ (∂A/∂p − λᵢ ∂B/∂p) rᵢ   with B-biorthogonal normalisation lᵢᴴ B rᵢ = 1
function TJLF._generalized_eigenvalues(A::Matrix{Complex{D}}, B::Matrix{Complex{D}}) where {D <: ForwardDiff.Dual}
    np  = ForwardDiff.npartials(D)
    Tag = ForwardDiff.tagtype(D)

    Af = map(a -> Complex{Float64}(ForwardDiff.value(real(a)), ForwardDiff.value(imag(a))), A)
    Bf = map(b -> Complex{Float64}(ForwardDiff.value(real(b)), ForwardDiff.value(imag(b))), B)

    # LAPACK ggev! returns left (L) and right (R) eigenvectors; it overwrites its
    # inputs, so pass copies.
    (alpha, beta, L, R) = ggev!('V', 'V', copy(Af), copy(Bf))
    λf = alpha ./ beta    # Vector{ComplexF64}
    n  = length(λf)

    # B-biorthogonal normalisation: rescale L so that lᵢᴴ B rᵢ = 1
    BfR = Bf * R
    for i in 1:n
        s = dot(L[:, i], BfR[:, i])
        if abs(s) > 1e-30
            L[:, i] ./= conj(s)
        end
    end

    dλ_re = zeros(n, np)
    dλ_im = zeros(n, np)
    for k in 1:np
        dAk = map(a -> Complex{Float64}(ForwardDiff.partials(real(a), k),
                                        ForwardDiff.partials(imag(a), k)), A)
        dBk = map(b -> Complex{Float64}(ForwardDiff.partials(real(b), k),
                                        ForwardDiff.partials(imag(b), k)), B)
        # Ck[i,i] = lᵢᴴ (dA − λᵢ dB) rᵢ
        LH_dAk_R = L' * (dAk * R)
        LH_dBk_R = L' * (dBk * R)
        for i in 1:n
            c = LH_dAk_R[i, i] - λf[i] * LH_dBk_R[i, i]
            dλ_re[i, k] = real(c)
            dλ_im[i, k] = imag(c)
        end
    end

    return map(1:n) do i
        Complex(ForwardDiff.Dual{Tag}(real(λf[i]), ntuple(k -> dλ_re[i, k], Val(np))...),
                ForwardDiff.Dual{Tag}(imag(λf[i]), ntuple(k -> dλ_im[i, k], Val(np))...))
    end
end

end
