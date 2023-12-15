# SpreadingDye

[![Build Status](https://github.com/mkborregaard/SpreadingDye.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mkborregaard/SpreadingDye.jl/actions/workflows/CI.yml?query=branch%3Amain)


This is a simple and fast implementation of the "spreading dye" (Jetz & Rahbek, 2001) algorithm for generating random geographical species ranges on a grid. The algorithm is implemented as conscientiously as possible from the description in the published article, without referencing any code, and as such is not guaranteed to behave exactly like the original analysis code.


Jetz, W., & Rahbek, C. (2001). Geometric constraints explain much of the species richness pattern in African birds. *Proceedings of the National Academy of Sciences*, **98**(10), 5661-5666.

### Exported functions:
```julia
spreading_dye(finalrange, dom, start = pick_random(dom); algo = rook)
```
Create a new random range with the spreading dye algorithm
- finalrange:   An integer specifying the number of grid cells to include in the range
- dom:          A boolean matrix or raster specifying the domain within which to put the range
- start:        The i, j position of the starting point of the range. If none specified it will take a random point on the domain
- algo:         The neighborhood that a cell can spread to. Either `rook`/ von Neumann neighborhood, or `queen` / Moore neighborhood

```julia
spreading_dye!(georange, finalrange, dom, start; algo = rook)
```
Create the random range in-place on an existing gird
- georange: the source matrix or raster file for the new range

```julia
center_cell(m::AbstractMatrix{<:Number})
```
Return the centroid cell of all cells with a non-0 finite value
- m:            A matrix or raster

```julia
pick_random(m)
```
Use rejection sampling to get a random point on the domain