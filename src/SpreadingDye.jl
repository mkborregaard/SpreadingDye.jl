module SpreadingDye

# the two types of neighborhood
const algos = Dict(
    :rook => ((-1,0), (0, 1), (1, 0), (0, -1)),
    :queen => Tuple((x,y) for x in -1:1, y in -1:1 if !(x == y == 0))    
)

function center_cell(m::AbstractMatrix{<:Number})
    range = 0
    x = 0
    y = 0
    for i in axes(m, 1), j in axes(m, 2)
        if isfinite(m[i,j]) && m[i, j] != 0
            x += i
            y += j
            range += 1
        end
    end
    round(Int, x/range), round(Int, y/range)
end

# is a point on the domain?
within_edges(dom::AbstractMatrix, point::Tuple{Int, Int}) = Base.checkbounds(Bool, dom, CartesianIndex(point))
on_domain(dom::AbstractMatrix{Bool}, point::Tuple{Int, Int}) = within_edges(dom, point) && dom[point...]

# use rejection sampling to get a random point on the domain
function random_point_on_domain(m::AbstractMatrix{Bool})
    x = (rand(axes(m,1)), rand(axes(m,2)))
    while !m[x...]
        x = (rand(axes(m,1)), rand(axes(m,2)))
    end
    x
end

# choose a random new cell
function _random_edge_or_jump(georange, edges, dom, algo)
    while !isempty(edges)
        newidx = rand(eachindex(edges))
        newcell = popat!(edges, newidx)
        !georange[newcell...] && return newcell
    end
    return jump(georange, dom, algo)
end

# grow the range one at a time
grow!(georange, edges, dom, algo::Symbol; ignore_domain = false) =
    _grow!(georange, edges, dom, algos[algo]; ignore_domain)

function _grow!(georange::AbstractMatrix{Bool}, edges::Vector, dom::AbstractMatrix{Bool}, algo::Tuple; ignore_domain = false)
    newcell = _random_edge_or_jump(georange, edges, dom, algo)
    for nb in algo
        neighbor = newcell .+ nb
        if within_edges(dom, neighbor) && 
            (ignore_domain || dom[neighbor...]) && 
            !georange[neighbor...]
            push!(edges, neighbor)
        end
    end
    georange[newcell...] = true
    newcell
end

# the internal range filling function
spreading_dye!(
    georange::AbstractMatrix{Bool}, 
    finalrange::Int, 
    dom::AbstractMatrix{Bool}, 
    start::Tuple{Int, Int} = random_point_on_domain(dom); 
    algo::Symbol = :rook
) = _spreading_dye!(georange, finalrange, dom, start, algos[algo])

function _spreading_dye!(
    georange::AbstractMatrix{Bool}, 
    finalrange::Int, 
    dom::AbstractMatrix{Bool}, 
    start::Tuple{Int, Int},
    algo::Tuple
)
    finalrange > count(dom) && error("finalrange $(finalrange) larger than domain size $(count(dom))")
    # is the start cell outside the domain?
    dom[start...] || return dom 

    fill!(georange, false) # TODO: Is this the right design? We might want to expand existing ranges
    # apply the starting cell
    georange[start...] = true
    edges = [(start .+ nb) for nb in algo if on_domain(dom, start .+ nb)]
    rangesize = 1

    # grow the range one cell at a time
    while (rangesize += 1) <= finalrange
         _grow!(georange, edges, dom, algo)
    end
    georange
end

# the outward_facing function to do the spreading dye
function spreading_dye(finalrange::Int, dom::AbstractMatrix{Bool}, start::Tuple{Int, Int} = random_point_on_domain(dom); algo::Symbol = :rook)
    georange = fill!(similar(dom, Bool), false)
    _spreading_dye!(georange, finalrange, dom, start, algos[algo])
end

function find_edges(georange::AbstractMatrix{Bool}, dom::AbstractMatrix{Bool}, algo::Tuple, ignore_domain = false)
    Tuple{Int,Int}[(i,j).+nb
    for i in axes(dom, 1), j in axes(dom, 2), nb in algo
        if within_edges(dom, (i,j).+nb) && georange[i,j] && !georange[((i,j).+nb)...] && (ignore_domain || dom[((i,j).+nb)...]) 
    ]
end

expand_spreading!(georange, add_cell, dom; algo::Symbol = :rook) =
    expand_spreading!(georange, add_cell, dom, algos[algo])
function expand_spreading!(georange::AbstractMatrix{Bool}, add_cells::Int, dom::AbstractMatrix{Bool}, algo::Tuple)
    add_cells + count(georange) > count(dom) && 
        error("not enough non-filled cells in domain to expand by $(add_cells) cells")
    edges = find_edges(georange, dom, algo)
    c = count(georange)
    for i in 1:add_cells
        _grow!(georange, edges, dom, algo)
        count(georange) == c + i || error("count(georange) should increase by $i, but is $(count(georange))")
    end
    georange
end

function jump(georange::AbstractMatrix{Bool}, dom::AbstractMatrix{Bool}, algo::Tuple)
    edges = find_edges(georange, dom, algo, true)
    newcell = rand(edges)
    while !dom[newcell...]
        newcell = _grow!(georange, edges, dom, algo; ignore_domain = true)
    end
    georange .&= dom
    georange[newcell...] = false
    return newcell
end

export spreading_dye, spreading_dye!, random_point_on_domain, center_cell, grow!, jump
end
