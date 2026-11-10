"""
Solver for the Pentomino Farm problem
    - Goal is to place the 12 pentominoes (https://en.wikipedia.org/wiki/Pentomino) on a square grid such that the maximal area is enclosed
    - One optimal example is given here: https://www.iread.it/lz/maximizing.html
    - Problem is converted into an ILP model (according to https://cs.stackexchange.com/questions/153818/maximize-enclosed-area-of-given-figures-on-2d-grid/153851#153851)
      and then solved with the commercial solver Gurobi, rectangular finite bounding box must be choosen
"""

using JuMP
using Gurobi
using ProgressMeter
using DelimitedFiles


"""
    solvePentomino(width, height, allSolutions)

Computes the optimal solution of the Pentomino Farm problem in the given bounding box,
the default 20x20 square is enough to cover all solutions (about ~350s). If allSolutions is true, all
5760 different solutions are explored (about ~750s).

# Arguments
- `width`: width of bounding box
- `height`: height of bounding box
"""
function solvePentomino(width=20, height=20; allSolutions=false)
    model = Model(Gurobi.Optimizer)

    # the 12 Pentominoes as binary arrays column by column (https://en.wikipedia.org/wiki/Pentomino)
    pieces = Vector{Matrix{Bool}}()
    push!(pieces, [false; true; false;; true; true; true;; true; false; false])  # F
    push!(pieces, [true;; true;; true;; true;; true])  # I
    push!(pieces, [true; true; true; true;; false; false; false; true])  # L
    push!(pieces, [false; true; true; true;; true; true; false; false])  # N
    push!(pieces, [true; true; true;; true; true; false])  # P
    push!(pieces, [true; false; false;; true; true; true;; true; false; false])  # T
    push!(pieces, [true; true;; false; true;; true; true])  # U
    push!(pieces, [true; true; true;; false; false; true;; false; false; true])  # V
    push!(pieces, [true; true; false;; false; true; true;; false; false; true])  # W
    push!(pieces, [false; true; false;; true; true; true;; false; true; false])  # X
    push!(pieces, [false; true; false; false;; true; true; true; true])  # Y
    push!(pieces, [true; false; false;; true; true; true;; false; false; true])  # Z
    
    @variable(model, x[1:width, 1:height], Bin)  # if tile is interior
    @variable(model, y[1:width, 1:height], Bin)  # if tile is boundary
    @variable(model, z[1:width, 1:height, 1:length(pieces), 1:8], Bin)  # z[w, h, p, d] the left top corner of pentomino p in orientation d is placed on (w, h) with (1, 1) being top left conrer

    @objective(model, Max, sum(x))  # Maximize area enclosed

    if allSolutions  # configure solver to find all solutions
        set_optimizer_attribute(model, "PoolSearchMode", 2)
        set_optimizer_attribute(model, "PoolSolutions", typemax(Int32))
        set_optimizer_attribute(model, "PoolGapAbs", 0.01)
    end

    #  Set of all z[w, h, p, d] that place a tile on (w, h)
    active = Dict{Pair{Int64, Int64}, Set{VariableRef}}()
    for i in 1 : width
        for j in 1 : height
            active[Pair(i, j)] = Set{VariableRef}()
        end
    end

    # every Pentomino is placed only once in one direction
    for i in 1 : length(pieces)
        @constraint(model, sum(z[:, :, i, :]) == 1)
    end

    # verify Pentominoes dont go out of bounds and fill active set
    for w in 1 : width
        for h in 1 : height
            for p in 1 : length(pieces)
                for d in 1 : 8
                    piece = rotate(pieces[p], d)

                    if (w + size(piece, 2) - 1 > width)  # too far right
                        @constraint(model, z[w, h, p, d] == 0)
                    elseif (h + size(piece, 1) - 1 > height)  # too far down
                        @constraint(model, z[w, h, p, d] == 0)
                    else
                        for i in 1 : size(piece, 2)
                            for j in 1 : size(piece, 1)
                                if (piece[j, i])
                                    push!(active[Pair(w + i - 1, h + j - 1)], z[w, h, p, d])
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    for w in 1 : width
        for h in 1 : height
            @constraint(model, sum(active[Pair(w, h)]) <= 1)  # no polyominoes overlap
            @constraint(model, sum(active[Pair(w, h)]) == y[w, h])  # border definition
        end
    end

    # if neighbour of interior is not boundary, than it is also interior: (not x[l]) or y[l'] or x[l']
    for w in 1 : width
        for h in 1 : height
            if (w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h] + x[w - 1, h] >= 1)  # left
            end
            if (h != 1)
                @constraint(model, (1 - x[w, h]) + y[w, h - 1] + x[w, h - 1] >= 1)  # above
            end
            if (w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h] + x[w + 1, h] >= 1)  # right
            end
            if (h != height)
                @constraint(model, (1 - x[w, h]) + y[w, h + 1] + x[w, h + 1] >= 1)  # below
            end

            if (w != 1 && h != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h - 1] + x[w - 1, h - 1] >= 1)  # left above
            end
            if (h != 1 && w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h - 1] + x[w + 1, h - 1] >= 1)  # right above
            end
            if (w != width && h != height)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h + 1] + x[w + 1, h + 1] >= 1)  # right below
            end
            if (h != height && w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h + 1] + x[w - 1, h + 1] >= 1)  # left below
            end
        end
    end
    
    # edges of box are neither boundary nor interior
    for w in 1 : width
        for h in 1 : height
            if (w == 1 || w == width || h == 1 || h == height)
                @constraint(model, x[w, h] == 0)
                @constraint(model, y[w, h] == 0)
            end
        end
    end

    # no tile is boundary and interior
    for w in 1 : width
        for h in 1 : height
            @constraint(model, x[w, h] + y[w, h] <= 1)
        end
    end

    # break translation symmetry by forcing fence in the top left corner
    @constraint(model, sum(y[2, :]) >= 1)
    @constraint(model, sum(y[:, 2]) >= 1)

    # break rotational symmetry by forcing F Pentomino in specific orientation
    @constraint(model, sum(z[:, :, 1, 1]) == 1)

    # exclude redundant orientations of rotational symmetric Pentominoes
    if allSolutions
        @constraint(model, sum(z[:, :, 2, 3:8]) == 0)  # I-Pentomino: 2 axis
        @constraint(model, sum(z[:, :, 6, 5:8]) == 0)  # T-Pentomino: 1 axis
        @constraint(model, sum(z[:, :, 7, 5:8]) == 0)  # U-Pentomino: 1 axis
        @constraint(model, sum(z[:, :, 8, 5:8]) == 0)  # V-Pentomino: 1 axis
        @constraint(model, sum(z[:, :, 9, 5:8]) == 0)  # W-Pentomino: 1 axis
        @constraint(model, sum(z[:, :, 10, 2:8]) == 0)  # X-Pentomino: 4 axis
        @constraint(model, sum(z[:, :, 12, 3:4]) == 0)  # Z-Pentomino: Rotational symmetry
        @constraint(model, sum(z[:, :, 12, 7:8]) == 0)
    end

    optimize!(model)

    if !allSolutions
        printSol(1, z, width, height, pieces)  # print first solution
        println(Int(round(objective_value(model))))
    else
        #for i in 1 : result_count(model)  # print all solutions
        #    printSol(i, z, width, height, pieces)
        #end
        println(Int(round(objective_value(model))))
        println()
        println("Number of solutions: " * string(result_count(model)))
        
        uniq = Set{Vector}()
        push!(uniq, order(1, z, width, height, pieces))

        printSol(1, z, width, height, pieces)

        output = hcat(Int64(sum(round.(value.(x)))), order(1, z, width, height, pieces)')
        @showprogress "Saving solutions to disk." for i in 1 : result_count(model)
            try
                ord = order(i, z, width, height, pieces)

                if !(ord in uniq)
                    output = [output; hcat(Int64(sum(round.(value.(x; result = i)))), ord')]
                    push!(uniq, ord)
                end
            catch y
                println(i)
                printSol(i, z, width, height, pieces)
            end
        end
        isdir("Results") || mkdir("Results")
        writedlm("Results/solutions-squ.csv",  output, ',')
    end
end


"""
    rotate(p, d)

For 1 <= d <= 4 the Pentomino p is rotated d-1 times counterclockwise. For 5 <= d <= 8 the Pentomino p is
rotated d-5 times counterclockwise and then mirrored left-right

# Arguments
- `p`: Pentomino
- `d`: direction
"""
function rotate(p::Matrix, d::Int64)
    if (d == 1)
        return p
    elseif (d == 2)
        return rotl90(p)
    elseif (d == 3)
        return rot180(p)
    elseif (d == 4)
        return rotr90(p)
    else
        return reverse(rotate(p, d - 4), dims=2)
    end  
end


"""
    printSol(n, z, width, height, pieces)

Print the n-th solution of the ILP model to the command line

# Arguments
- `n`: number of solution
- `z`: variable matrix of z
- `width`: width of bounding box
- `height`: height of bounding box
- `pieces`: Vector of Pentominoes
"""
function printSol(n::Int64, z, width, height, pieces)  # print n-th solution
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))
    for i in 1 : length(solution)
        pieceSet = rotate(pieces[solution[i][3]], solution[i][4])

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end

    display(output)
end


"""
    order(n, z, width, height, pieces)

Find the order of the pieces in the n-th solution. For uniqueness the first Pentomino in the output is always 1
and the second piece has a bigger number then the last piece, so for example [1, 3, 10, 5, 8, 7, 4, 12, 11, 6, 9, 2]
and not [1, 2, 9, 6, 11, 12, 4, 7, 8, 5, 10, 3].

# Arguments
- `n`: number of solution
- `z`: variable matrix of z
- `width`: width of bounding box
- `height`: height of bounding box
- `pieces`: Vector of Pentominoes
"""
function order(n::Int64, z, width, height, pieces)  # determine order of pieces
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))  # create output matrix
    for i in 1 : length(solution)
        pieceSet = rotate(pieces[solution[i][3]], solution[i][4])

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end

    neighList = Vector{Set{Int}}()  # find two neighbours of every pentomino
    for i in 1 : length(pieces)
        push!(neighList, Set{Int}())
    end

    for x in 1 : width
        for y in 1 : height
            for d in ((0,1),(0,-1),(1,0),(-1,0))
                neigh = (x, y) .+ d

                if (neigh[1] > 0 && neigh[2] > 0 && neigh[1] <= width && neigh[2] <= height)
                    if output[x, y] > 0 && output[neigh[1], neigh[2]] > 0 && output[x, y] != output[neigh[1], neigh[2]]
                        push!(neighList[output[x, y]], output[neigh[1], neigh[2]])
                    end
                end
            end
        end
    end

    ord = Vector{Int}()  # bring neighbour list in order
    push!(ord, 1)
    discovered = Set{Int}()
    push!(discovered, 1)
    while !isempty(setdiff(neighList[ord[length(ord)]], discovered))
        n = maximum(setdiff(neighList[ord[length(ord)]], discovered))
        push!(discovered, n)
        push!(ord, n)
    end

    return ord
end
