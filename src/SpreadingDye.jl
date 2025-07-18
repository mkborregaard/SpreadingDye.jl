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
within_edges(dom::AbstractMatrix, point::Tuple{Int, Int}) = min(point...) > 0 && first(point) <= size(dom, 1) && last(point) <= size(dom, 2)
on_domain(dom::AbstractMatrix{Bool}, point::Tuple{Int, Int}) = within_edges(dom, point) && dom[point...]
# is a point an edge of the domain?
is_edge(georange, dom, point, ignore_domain) = 
    within_edges(dom, point) && (ignore_domain || dom[point...]) && !georange[point...]

# use rejection sampling to get a random point on the domain
function random_point_on_domain(m::AbstractMatrix{Bool})
    x = (rand(axes(m,1)), rand(axes(m,2)))
    while !m[x...]
        x = (rand(axes(m,1)), rand(axes(m,2)))
    end
    x
end

# grow the range one at a time
function grow!(georange::AbstractMatrix{Bool}, edges::Set, dom::AbstractMatrix{Bool}, algo::Tuple; ignore_domain = false)
    newcell = isempty(edges) ? random_point_on_domain(georange) : first(rand(edges))
    while georange[newcell...]
        isempty(edges) && push!(edges, jump(georange, dom, algo)) # allows for patchy ranges
        edge, newcell = rand(edges)
        pop!(edges, (edge, newcell))
    end
    georange[newcell...] = true
    for nb in algo
        neighbor = newcell .+ nb
        if is_edge(georange, dom, neighbor, ignore_domain)
            push!(edges, (newcell, neighbor))
        end
    end
    newcell
end

# the internal range filling function
spreading_dye!(georange, finalrange, dom, start = random_point_on_domain(dom); algo::Symbol = :rook) =
    spreading_dye!(georange, finalrange, dom, start, algos[algo])
function spreading_dye!(georange::AbstractMatrix{Bool}, finalrange::Int, dom::AbstractMatrix{Bool}, start::Tuple{Int, Int}, algo::Tuple)
    finalrange > count(dom) && error("finalrange $(finalrange) larger than domain size $(count(dom))")
    # is the start cell outside the domain?
    dom[start...] || return dom 

    fill!(georange, false) # TODO: Is this the right design? We might want to expand existing ranges
    # apply the starting cell
    georange[start...] = true
    edges = Set((start, start .+ nb) for nb in algo if on_domain(dom, start .+ nb))
    rangesize = 1

    # grow the range one cell at a time
    while (rangesize += 1) <= finalrange
         grow!(georange, edges, dom, algo)
    end
    georange
end

# the outward_facing function to do the spreading dye
function spreading_dye(finalrange::Int, dom::AbstractMatrix{Bool}, start::Tuple{Int, Int} = random_point_on_domain(dom); algo::Symbol = :rook)
    georange = fill!(similar(dom, Bool), false)
    spreading_dye!(georange, finalrange, dom, start; algo)
end

function find_edges(georange::AbstractMatrix{Bool}, dom::AbstractMatrix{Bool}, algo, ignore_domain::Bool = false)
    Set(((i,j), (i,j).+nb) 
    for i in axes(dom, 1), j in axes(dom, 2), nb in algo
        if georange[i,j] && is_edge(georange, dom, (i,j).+nb, ignore_domain)
    )
end

expand_spreading!(georange, add_cells, dom; algo::Symbol = :rook) = 
    expand_spreading!(georange, add_cells, dom, algos[algo])
function expand_spreading!(georange::AbstractMatrix{Bool}, add_cells::Int, dom::AbstractMatrix{Bool}, algo::Tuple)
    add_cells + count(georange) > count(dom) && 
        error("not enough non-filled cells in domain to expand by $(add_cells) cells")
    edges = find_edges(georange, dom, algo)
    for i in 1:add_cells
        grow!(georange, edges, dom, algo)
    end
    georange
end

const i = [0]
function jump(georange::AbstractMatrix{Bool}, dom::AbstractMatrix{Bool}, algo)
    edges = find_edges(georange, dom, algo, true)
    newcell = last(rand(edges))
    while !dom[newcell...]
        newcell = grow!(georange, edges, dom, algo; ignore_domain = true)
    end
    georange .&= dom
    georange[newcell...] = false
    ((0,0), newcell)
end

export spreading_dye, spreading_dye!, random_point_on_domain, center_cell, grow!, jump
end
