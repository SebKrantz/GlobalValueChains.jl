# Core data type and precomputation for ICIO decompositions.
#
# Layout convention (matches the icio / decompr ordering): country-major, sector-minor.
# Country-sector index for country g (1..G), sector n (1..N) is (g-1)*N + n.

"""
    ICIOModel

Holds an Inter-Country Input-Output table and all matrices/vectors precomputed once and
reused across decompositions. Construct with [`load_icio`](@ref) or [`read_icio_csv`](@ref).

Fields of interest: `G` (countries), `N` (sectors), `GN = G*N`, `regions`, `sectors`,
`X` (output), `V` (value-added coefficients), `A`, `B` (global Leontief inverse), `L`
(block-diagonal local Leontief), `E` (total exports), `ESR` (bilateral exports by destination,
`GN×G`), `FD` (final demand by absorbing country, `GN×G`).
"""
struct ICIOModel
    G::Int
    N::Int
    GN::Int
    regions::Vector{String}
    sectors::Vector{String}
    X::Vector{Float64}        # gross output (GN)
    V::Vector{Float64}        # value-added coefficients v/o (GN)
    A::Matrix{Float64}        # input coefficients T ./ X' (GN×GN)
    B::Matrix{Float64}        # global Leontief inverse (I-A)^{-1} (GN×GN)
    L::Matrix{Float64}        # block-diagonal local Leontief (I-A_ss)^{-1} per country
    E::Vector{Float64}        # total exports per country-sector (GN)
    ESR::Matrix{Float64}      # bilateral exports by destination country (GN×G)
    FD::Matrix{Float64}       # final demand by absorbing country = Y (GN×G)
    # cached coefficients reused by the source/exporter decompositions
    VBdom::Vector{Float64}    # Σ_{i∈ctry(j)} V_i B_ij           (domestic VA multiplier)
    VBfor::Vector{Float64}    # Σ_{i∉ctry(j)} V_i B_ij           (foreign  VA multiplier)
    VLdom::Vector{Float64}    # Σ_{i∈ctry(j)} V_i L_ij           (local domestic VA multiplier)
    fvacoef::Vector{Float64}  # exporter/source foreign-VA-once coefficient (GN)
end

# country index of country-sector i (1-based)
@inline ctry(i::Integer, N::Integer) = (i - 1) ÷ N + 1
# 1:GN index range of country g's block
@inline blockrange(g::Integer, N::Integer) = (g - 1) * N + 1 : g * N

"""
    load_icio(VA, FD, T; regions, sectors, X = nothing)

Build an [`ICIOModel`](@ref) from the three core matrices, mirroring `decompr`'s
`load_tables_vectors`:

* `T`  — `GN×GN` intermediate transactions (`T[i,j]` = inputs from country-sector `i` used by `j`).
* `FD` — `GN×G` final demand (`FD[i,r]` = final goods `i` absorbed in country `r`).
* `VA` — length-`GN` value added, or `nothing` to use the icio residual `X .- vec(sum(T;dims=1))`.

`regions` (length `G`) and `sectors` (length `N`) are country/industry names. `X` (output) is
computed as `rowSums(T) + rowSums(FD)` when not supplied.

Note: to reproduce `icio` exactly, value added is the column residual of the table. Passing a
`VA` that differs from `X .- colSums(T)` makes column sums of `V*B` deviate from 1 (so
`DC + FC` may differ slightly from gross exports) — faithful to the supplied data.
"""
function load_icio(VA, FD::AbstractMatrix, T::AbstractMatrix;
                   regions::AbstractVector, sectors::AbstractVector,
                   X::Union{Nothing,AbstractVector} = nothing)
    G = length(regions)
    N = length(sectors)
    GN = G * N
    size(T) == (GN, GN) || throw(DimensionMismatch(
        "T is $(size(T)), expected ($GN, $GN) = (G*N, G*N)"))
    size(FD, 1) == GN || throw(DimensionMismatch(
        "FD has $(size(FD,1)) rows, expected $GN"))
    size(FD, 2) == G || throw(DimensionMismatch(
        "FD has $(size(FD,2)) columns, expected G = $G (one final-demand column per country)"))

    Tm = Matrix{Float64}(T)
    FDm = Matrix{Float64}(FD)

    Xv = X === nothing ? vec(sum(Tm; dims = 2)) .+ vec(sum(FDm; dims = 2)) :
                         Vector{Float64}(X)
    length(Xv) == GN || throw(DimensionMismatch("X has length $(length(Xv)), expected $GN"))

    # value added: supplied, else icio residual (= X - colSums(T))
    Vav = VA === nothing ? Xv .- vec(sum(Tm; dims = 1)) : Vector{Float64}(VA)
    length(Vav) == GN || throw(DimensionMismatch("VA has length $(length(Vav)), expected $GN"))

    # A = T column-normalised by output; guard X==0
    A = similar(Tm)
    @inbounds for j in 1:GN
        xj = Xv[j]
        invx = xj == 0 ? 0.0 : 1.0 / xj
        for i in 1:GN
            A[i, j] = Tm[i, j] * invx
        end
    end

    Vc = similar(Xv)
    @inbounds for i in 1:GN
        Vc[i] = Xv[i] == 0 ? 0.0 : Vav[i] / Xv[i]
    end

    # Global Leontief inverse B = (I - A)^{-1}  -- the single GN×GN inversion
    ImA = Matrix{Float64}(I, GN, GN)
    ImA .-= A
    B = ImA \ Matrix{Float64}(I, GN, GN)

    # Block-diagonal local Leontief L: blocks (I - A_ss)^{-1}
    L = zeros(Float64, GN, GN)
    for g in 1:G
        rng = blockrange(g, N)
        Lss = (Matrix{Float64}(I, N, N) .- @view A[rng, rng]) \ Matrix{Float64}(I, N, N)
        @view(L[rng, rng]) .= Lss
    end

    # Exports: total E (GN) and bilateral ESR (GN×G), zeroing within-country use
    E = zeros(Float64, GN)
    ESR = zeros(Float64, GN, G)
    @inbounds for g in 1:G
        rrows = blockrange(g, N)
        for i in rrows
            tot = 0.0
            for r in 1:G
                r == g && continue
                crng = blockrange(r, N)
                s = 0.0
                for j in crng
                    s += Tm[i, j]
                end
                s += FDm[i, r]
                ESR[i, r] = s
                tot += s
            end
            E[i] = tot
        end
    end

    # cached coefficients: VBdom/VBfor/VLdom and exporter/source fvacoef
    VBdom = zeros(Float64, GN)
    VBfor = zeros(Float64, GN)
    VLdom = zeros(Float64, GN)
    fvacoef = zeros(Float64, GN)
    II_N = Matrix{Float64}(I, N, N)
    for g in 1:G
        rng = blockrange(g, N)
        # column sums of V.*B over all rows (for VBfor) and over domestic block (VBdom)
        @inbounds for j in rng
            dom = 0.0
            for i in rng
                dom += Vc[i] * B[i, j]
            end
            total = 0.0
            for i in 1:GN
                total += Vc[i] * B[i, j]
            end
            VBdom[j] = dom
            VBfor[j] = total - dom
            ld = 0.0
            for i in rng
                ld += Vc[i] * L[i, j]
            end
            VLdom[j] = ld
        end
        # M_g = Σ_{j≠g} A_gj B_jg = (A[g,:] * B[:,g]) - A_gg B_gg     (N×N)
        Ag = @view A[rng, :]
        Bg = @view B[:, rng]
        Mg = Ag * Bg
        Mg .-= @view(A[rng, rng]) * @view(B[rng, rng])
        # fvacoef[g-block] = (I+M_g)^{-T} VBfor[g-block]
        IM = II_N .+ Mg
        fvacoef[rng] .= transpose(IM) \ VBfor[rng]
    end

    return ICIOModel(G, N, GN, collect(String, regions), collect(String, sectors),
                     Xv, Vc, A, B, L, E, ESR, FDm, VBdom, VBfor, VLdom, fvacoef)
end

function Base.show(io::IO, m::ICIOModel)
    print(io, "ICIOModel(", m.G, " countries × ", m.N, " sectors = ", m.GN,
          " country-sectors)")
end
