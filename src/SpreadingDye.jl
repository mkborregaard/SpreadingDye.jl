module SpreadingDye

# the two types of neighborhood
const rook = ((-1,0), (0, 1), (1, 0), (0, -1))
const quen = Tuple((x,y) for x in -1:1, y in -1:1 if !(x == y == 0))

# the internal range filling function
function spreading_dye!(georange, finalrange, dom, start; algo = rook)
    # is the start cell outside the domain?
    dom[start...] || return dom 

    fill!(georange, false)
    # apply the starting cell
    georange[start...] = true
    edges = Set((start, start .+ nb) for nb in algo if on_domain(dom, start .+ nb))
    rangesize = 1
    newcell = edge = start

    # grow the range one cell at a time
    while (rangesize += 1) <= finalrange
        grow!(georange, edges, dom, algo)
    end
end

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
on_domain(dom, point) = min(point...) > 0 && first(point) <= size(dom, 1) && last(point) <= size(dom, 2) && dom[point...]

# use rejection sampling to get a random point on the domain
function pick_random(dom)
    x = (rand(axes(dom,1)), rand(axes(dom,2)))
    while !dom[x...]
        x = (rand(axes(dom,1)), rand(axes(dom,2)))
    end
    x
end

# the outward_facing function to do the spreading dye
function spreading_dye(finalrange, dom, start = pick_random(dom); algo = rook)
    georange = fill!(similar(dom, Bool), false)
    spreading_dye!(georange, finalrange, dom, start; algo = rook)
    georange
end

# grow the range one at a time
function grow!(georange, edges, dom, algo)
    edge, newcell = rand(edges)
    pop!(edges, (edge, newcell))
    georange[newcell...] = true
    for nb in algo
        neighbor = newcell .+ nb
        on_domain(dom, neighbor) && !georange[neighbor...] && push!(edges, (newcell, neighbor))
    end
end

export spreading_dye, spreading_dye!, pick_random, center_cell
end
