using ICIO
using DataFrames
using Test
using LinearAlgebra: I

# -------------------------------------------------------------------------
# Build a small but economically valid synthetic ICIO table (deterministic).
# G countries, N sectors. Intermediate use kept small so value added > 0.
# -------------------------------------------------------------------------
function toy_table(; G = 3, N = 2, seed = 1)
    GN = G * N
    # simple deterministic positive entries
    T = [0.4 + 0.01 * ((i * 7 + j * 3 + seed) % 11) for i in 1:GN, j in 1:GN]
    FD = [5.0 + 0.1 * ((i * 5 + r * 2 + seed) % 13) for i in 1:GN, r in 1:G]
    regions = ["C$g" for g in 1:G]
    sectors = ["S$n" for n in 1:N]
    return T, FD, regions, sectors
end

approxcols(df, a, b; atol = 1e-7) = maximum(abs.(Float64.(df[!, a]) .- Float64.(df[!, b]))) < atol

@testset "ICIO.jl" begin
    T, FD, regions, sectors = toy_table(G = 4, N = 3)
    G, N = length(regions), length(sectors)

    @testset "construction" begin
        m = load_icio(nothing, FD, T; regions = regions, sectors = sectors)
        @test m.G == G && m.N == N && m.GN == G * N
        @test size(m.B) == (G * N, G * N)
        # residual VA ⇒ column sums of V*B == 1 ⇒ DC + FC == gross exports
        @test all(c -> abs(c - 1) < 1e-9, vec(sum(m.V .* m.B; dims = 1)))
        @test_throws DimensionMismatch load_icio(nothing, FD, T;
            regions = regions[1:end-1], sectors = sectors)
    end

    m = load_icio(nothing, FD, T; regions = regions, sectors = sectors)

    @testset "sector identities (exporter/source)" begin
        s = decompose(m; level = :sector)
        @test nrow(s) == G * N
        @test approxcols(s, :gexp, :dc) == false   # sanity: not trivially equal
        for (lhs, parts) in (
                (:gexp, [:dc, :fc]), (:dc, [:dva, :ddc]), (:fc, [:fva, :fdc]),
                (:dva, [:vax, :ref]), (:gvc, [:gvcb, :gvcf]))
            @test maximum(abs.(s[!, lhs] .- reduce(+, (s[!, p] for p in parts)))) < 1e-9
        end
        @test maximum(abs.(s.gvc .- (s.gexp .- s.davax))) < 1e-9
        @test maximum(abs.(s.gvcb .- (s.fc .+ s.ddc))) < 1e-9
    end

    @testset "additivity bilateral → sector → country" begin
        sec = decompose(m; level = :sector)
        bil = decompose(m; level = :bilateral)
        @test nrow(bil) == G * (G - 1) * N
        bagg = combine(groupby(bil, [:from_region, :from_sector]),
                       [:gexp, :dva, :fva, :gvc] .=> sum .=> [:gexp, :dva, :fva, :gvc])
        sec_k = sort(sec, [:from_region, :from_sector])
        bil_k = sort(bagg, [:from_region, :from_sector])
        for c in [:gexp, :dva, :fva, :gvc]
            @test maximum(abs.(sec_k[!, c] .- bil_k[!, c])) < 1e-8
        end
        # country source == aggregate(sector)
        cty = decompose(m; level = :country, perspective = :exporter, approach = :source)
        cagg = combine(groupby(transform(sec, :from_region), :from_region),
                       :gexp => sum => :gexp)
        @test maximum(abs.(sort(cty, :country).gexp .- sort(cagg, :from_region).gexp)) < 1e-8
    end

    @testset "world/sink vs exporter/source share domestic side + FC" begin
        w = decompose(m; level = :country, perspective = :world, approach = :sink)
        s = decompose(m; level = :country, perspective = :exporter, approach = :source)
        for c in [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc]
            @test maximum(abs.(w[!, c] .- s[!, c])) < 1e-8
        end
        # both split FC into FVA + FDC
        @test maximum(abs.(w.fc .- (w.fva .+ w.fdc))) < 1e-8
        @test ncol(w) == 1 + 9    # country + 9 terms
        @test ncol(s) == 1 + 13   # country + 13 terms
    end

    @testset "read_icio_csv round-trip" begin
        dir = mktempdir()
        tpath = joinpath(dir, "table.csv")
        cpath = joinpath(dir, "clist.csv")
        open(tpath, "w") do io
            M = hcat(T, FD)
            for i in 1:size(M, 1)
                println(io, join(M[i, :], ","))
            end
        end
        open(cpath, "w") do io
            for r in regions; println(io, r); end
        end
        mc = read_icio_csv(tpath, cpath; sectors = sectors)
        @test mc.G == G && mc.N == N
        @test maximum(abs.(mc.B .- m.B)) < 1e-10
        dc1 = decompose(mc; level = :country, perspective = :world, approach = :sink)
        dc2 = decompose(m; level = :country, perspective = :world, approach = :sink)
        @test maximum(abs.(dc1.fva .- dc2.fva)) < 1e-9
    end

    @testset "multi-year batch" begin
        m2 = load_icio(nothing, FD .* 1.1, T .* 0.9; regions = regions, sectors = sectors)
        d = decompose(Dict(2020 => m, 2021 => m2); level = :country)
        @test "year" in names(d)
        @test sort(unique(d.year)) == [2020, 2021]
        @test nrow(d) == 2 * G
    end

    @testset "option validation" begin
        @test_throws ErrorException decompose(m; level = :sector,
            perspective = :world, approach = :sink)
        @test_throws ErrorException decompose(m; level = :nonsense)
    end
end
