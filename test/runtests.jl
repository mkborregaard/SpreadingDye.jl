using SpreadingDye
using Test

@testset "Simple example" begin
    # create a new domain and fill it
    dom = fill(true, 100, 100)
    start = (50,50)
    sd = spreading_dye(3000, dom, start)

end

@testset "Jesper's example" begin
    using Rasters
    dd = Raster("../data/domaene.tif")
    dom = !=(missingval(dd)).(dd)
    rr = Raster("../data/species_range.tif")
    rangesize = count(isfinite, rr)
    coords = collect(CartesianIndices(rr))[isfinite.(rr)]
    centroid = center_cell(rr)
    sd = spreading_dye(rangesize, dom, centroid)
end

@testset "South America birds empirical data" begin
    # load libraries
    using CSV, DataFrames
    bd = CSV.read("../data/Samerica_birds_16Apr2010.txt", DataFrame, header = false)

    # reshape data nand initialize objects
    xa = extrema(bd.Column4)
    ya = extrema(bd.Column5)
    d = Dict{String, Int}()
    longs = Dict{String, Int}()
    lats = Dict{String, Int}()
    dom = falses(Int(last(xa) - first(xa) + 2), Int(last(ya) - first(ya) + 2))
    for row in eachrow(bd)
        x = Int(row.Column4 - first(xa) + 2)
        y = Int(row.Column5 - first(ya) + 2)
        dom[x, y] = true
        if !haskey(d, row.Column1)
            d[row.Column1] = 0
            longs[row.Column1] = 0
            lats[row.Column1] = 0
        end
        d[row.Column1] += 1
        longs[row.Column1] += x
        lats[row.Column1] += y
    end

    centroids = Dict(sp => (longs[sp] ÷ d[sp], lats[sp] ÷  d[sp]) for sp in keys(d))

    # do a spreading dye for the first species
    sp = first(keys(d))
    sd = spreading_dye(d[sp], dom, centroids[sp])
    heatmap(range(xa...), range(ya...), sd)

    # do it for all
    function doall(dom, d, centroids)
        rich = zeros(size(dom))
        georange = fill(false, size(dom))
        for sp in keys(d)
            spreading_dye!(georange, d[sp], dom, centroids[sp])
            rich .+= georange
        end
        rich[.!(dom)] .= NaN
        rich
    end

    doall(dom, d, centroids)
end