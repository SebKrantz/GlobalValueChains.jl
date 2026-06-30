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

    # --- helpers reused below ---
    identities(df) = all((
        maximum(abs.(df.gexp .- (df.dc .+ df.fc))) < 1e-8,
        maximum(abs.(df.dc .- (df.dva .+ df.ddc))) < 1e-8,
        maximum(abs.(df.fc .- (df.fva .+ df.fdc))) < 1e-8))
    nonneg(df, cols) = all(minimum(df[!, c]) > -1e-8 for c in cols)
    byctry(df, cols) = sort(combine(groupby(df, :from_region),
                                    (cols .=> sum .=> cols)...), :from_region)

    @testset "world source/sink (country, 9 terms)" begin
        wsk = decompose(m; perspective = :world, approach = :sink)
        wso = decompose(m; perspective = :world, approach = :source)
        src = decompose(m; perspective = :exporter, approach = :source)
        for w in (wsk, wso)
            @test ncol(w) == 1 + 9
            @test identities(w)
            @test nonneg(w, [:dva, :vax, :ref, :ddc, :fva, :fdc])
            # domestic side + FC are perspective-invariant (== exporter/source)
            for c in [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc]
                @test maximum(abs.(w[!, c] .- src[!, c])) < 1e-8
            end
        end
        # world source & sink FVA differ by country but match in the world total (BM2019 §5.1)
        @test abs(sum(wso.fva) - sum(wsk.fva)) < 1e-7
        @test maximum(abs.(wso.fva .- wsk.fva)) > 1e-6      # genuinely different per country
    end

    @testset "exporter/sink sector & bilateral (FVA uniqueness, eq. 33–39)" begin
        sks = decompose(m; level = :sector, approach = :sink)
        sss = decompose(m; level = :sector, approach = :source)
        skb = decompose(m; level = :bilateral, approach = :sink)
        ssb = decompose(m; level = :bilateral, approach = :source)
        sc  = sort(decompose(m; level = :country), :country)
        @test ncol(sks) == 2 + 9
        @test ncol(skb) == 3 + 10           # adds vaxim
        @test identities(sks) && identities(skb)
        @test nonneg(sks, [:dva, :vax, :ref, :ddc, :fva, :fdc])
        @test nonneg(skb, [:dva, :vax, :vaxim, :ref, :ddc, :fva, :fdc])
        # DC and FC are identical to the source breakdown at every cell (BM2019 §3.2)
        @test maximum(abs.(sks.dc .- sss.dc)) < 1e-10
        @test maximum(abs.(skb.fc .- ssb.fc)) < 1e-10
        # anchor: summed over importers (& sectors) sink == source country (Σ_r FVAsink = Σ_r FVAsource)
        aggs = byctry(sks, [:dva, :fva, :vax, :ref])
        for c in [:dva, :fva, :vax, :ref]
            @test maximum(abs.(aggs[!, c] .- sc[!, c])) < 1e-7
        end
        # bilateral sink: additivity to sector, and VAXIM nests DAVAX ⊆ VAXIM ⊆ VAX
        aggb = byctry(skb, [:dva, :fva])
        @test maximum(abs.(aggb.dva .- sc.dva)) < 1e-7
        @test minimum(skb.vaxim .- ssb.davax) > -1e-9    # DAVAX ⊆ VAXIM
        @test minimum(skb.vax .- skb.vaxim) > -1e-9      # VAXIM ⊆ VAX
    end

    @testset "self perimeter (sectexp / sectbil, 7 terms)" begin
        se = decompose(m; level = :sector, perspective = :self)
        sb = decompose(m; level = :bilateral, perspective = :self)
        sss = decompose(m; level = :sector, approach = :source)
        ssb = decompose(m; level = :bilateral, approach = :source)
        skb = decompose(m; level = :bilateral, approach = :sink)
        @test ncol(se) == 2 + 7 && ncol(sb) == 3 + 7
        @test identities(se) && identities(sb)
        @test nonneg(se, [:dva, :ddc, :fva, :fdc]) && nonneg(sb, [:dva, :ddc, :fva, :fdc])
        # content is perimeter-invariant; DVA★/FVA★ dominate the country-perimeter measures (eq. 46)
        @test maximum(abs.(se.dc .- sss.dc)) < 1e-10
        @test minimum(sb.dva .- ssb.dva) > -1e-9
        @test minimum(sb.dva .- skb.dva) > -1e-9
        @test minimum(sb.fva .- ssb.fva) > -1e-9
    end

    @testset "imports (importer perspective, eq. 51)" begin
        ic = decompose(m; flow = :imports)
        ib = decompose(m; flow = :imports, level = :bilateral)
        @test names(ic) == ["importer", "gimp", "va", "dc"]
        @test names(ib) == ["importer", "origin", "va", "dc"]
        @test maximum(abs.(ic.gimp .- (ic.va .+ ic.dc))) < 1e-8
        @test all(ic.va .> -1e-8) && all(ic.dc .> -1e-8)
        @test all(ib.va .> -1e-8) && all(ib.dc .> -1e-8)
        # gross imports + world consistency: Σ_r imports == Σ exports
        totexp = sum(decompose(m; level = :country).gexp)
        @test abs(sum(ic.gimp) - totexp) < 1e-7
        # by-origin VA/DC sum (over origins) to the country-level totals
        aggi = sort(combine(groupby(ib, :importer), [:va, :dc] .=> sum .=> [:va, :dc]), :importer)
        ics = sort(ic, :importer)
        @test maximum(abs.(aggi.va .- ics.va)) < 1e-7
        @test maximum(abs.(aggi.dc .- ics.dc)) < 1e-7
    end

    @testset "option validation" begin
        @test_throws ErrorException decompose(m; level = :sector, perspective = :world)
        @test_throws ErrorException decompose(m; level = :bilateral, perspective = :world)
        @test_throws ErrorException decompose(m; level = :country, perspective = :self)
        @test_throws ErrorException decompose(m; flow = :imports, level = :sector)
        @test_throws ErrorException decompose(m; flow = :nonsense)
        @test_throws ErrorException decompose(m; level = :nonsense)
        @test_throws ErrorException decompose(m; level = :sector, approach = :nonsense)
    end

    @testset "multi-year batch with flow/approach" begin
        m2 = load_icio(nothing, FD .* 1.1, T .* 0.9; regions = regions, sectors = sectors)
        d = decompose(Dict(2020 => m, 2021 => m2); level = :bilateral, approach = :sink)
        @test "year" in names(d) && "vaxim" in names(d)
        @test nrow(d) == 2 * G * (G - 1) * N
        di = decompose(Dict(2020 => m, 2021 => m2); flow = :imports)
        @test sort(unique(di.year)) == [2020, 2021] && "gimp" in names(di)
    end
end
