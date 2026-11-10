using JuMP
using Gurobi
using DelimitedFiles
using ProgressMeter

function prettyPrint(p::Matrix{Int64})
    for i in 1 : size(p, 1)
        if i % 2 == 1
            print("[")
        else
            print(" [")
        end

        for j in 1 : size(p, 2)
            if p[i, j] != 0
                print(Char(p[i, j] - 1 + 'A'))
            else
                print(" ")
            end
            print("][")
        end
        println("")
    end
end

function prettyPrint(p::Matrix{Bool})
    for i in 1 : size(p, 1)
        if i % 2 == 1
            print("[")
        else
            print(" [")
        end

        for j in 1 : size(p, 2)
            if p[i, j]
                print("*")
            else
                print(" ")
            end
            print("][")
        end
        println("")
    end
end

function printSol(n::Int64, z, width, height, pieces)  # print n-th solution
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))
    for i in 1 : length(solution)
        pieceSet = pieces[solution[i][3]][solution[i][4]]

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end

    prettyPrint(output)
end

function order(n::Int64, z, width, height, pieces)  # determine order of pieces
    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z; result = n)))  # create output matrix
    for i in 1 : length(solution)
        pieceSet =  pieces[solution[i][3]][solution[i][4]]

        for j in 1 : size(pieceSet, 1)
            for l in 1 : size(pieceSet, 2)
                if (pieceSet[j, l])
                    output[j + solution[i][2] - 1, l + solution[i][1] - 1] = solution[i][3]
                end
            end
        end
    end
    
    neighList = Vector{Set{Int}}()  # find two neighbours of every polyhex
    for i in 1 : length(pieces)
        push!(neighList, Set{Int}())
    end

    for x in 1 : width
        for y in 1 : height
            for d in ((-1, 0), (0, -1), (1, 0), (0, 1), (1, 1), (1, -1), (-1, 1), (-1, -1))
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

function solvePolyhex(area = Pair(16, 16))  # TEMP 16
    """
    Setup
    """

    # bounding box to search in
    width = area[1]
    height = area[2]

    pieces = Vector{Vector{Matrix{Bool}}}()
    for i in 1 : 7
        push!(pieces, Vector{Matrix{Bool}}())
    end

    # petrahexes (https://en.wikipedia.org/wiki/Polyhex_(mathematics))
    push!(pieces[1], [true;; true;; true;; true])  # bar
    push!(pieces[1], [true; true; false; false;; false; false; true; true])
    push!(pieces[1], [false; false; false; true;; false; true; true; false;; true; false; false; false])

    push!(pieces[2], [false; true;; true; false;; true; false;; true; false])  # worm
    push!(pieces[2], [false; true;; false; true;; false; true;; true; false])
    push!(pieces[2], [true; true;; false; true;; false; true])
    push!(pieces[2], [true; false;; true; false;; true; true])
    push!(pieces[2], [true; true; false; true;; false; false; true; false])
    push!(pieces[2], [false; true; false; false;; false; true; false; false;; false; false; true; true])
    push!(pieces[2], [true; true; false;; false; false; true;; false; false; true])
    push!(pieces[2], [false; true; false; false;; true; false; true; true])
    push!(pieces[2], [false; false; false; true;; true; true; true; false])
    push!(pieces[2], [false; true; true; true;; true; false; false; false])
    push!(pieces[2], [false; true; true;; true; false; false;; true; false; false])
    push!(pieces[2], [false; false; true; false;; false; true; true; false;; true; false; false; false])

    push!(pieces[3], [true; true;; true; false;; true; false])  # pistol
    push!(pieces[3], [true; false;; true; true;; true; false])
    push!(pieces[3], [false; true;; true; true;; false; true])
    push!(pieces[3], [false; true;; false; true;; true; true])
    push!(pieces[3], [true; true; true;; true; false; false])
    push!(pieces[3], [false; true; true;; true; true; false])
    push!(pieces[3], [false; true; false;; false; true; true;; true; false; false])
    push!(pieces[3], [false; true; true;; true; false; true])
    push!(pieces[3], [false; true; false;; true; true; false;; false; false; true])
    push!(pieces[3], [true; true; true;; false; false; true])
    push!(pieces[3], [true; true; false;; false; true; true])
    push!(pieces[3], [true; true; false;; true; false; true])
    
    push!(pieces[4], [true; true; true;; false; true; false])  # propellor

    push!(pieces[5], [true; true;; false; true;; true; false])  # arch
    push!(pieces[5], [false; true;; true; false;; true; true])
    push!(pieces[5], [false; false; true;; true; true; true])
    push!(pieces[5], [false; true; false;; true; false; true;; true; false; false])
    push!(pieces[5], [true; false; false;; true; true; true])
    push!(pieces[5], [false; true; false;; true; false; true;; false; false; true])

    push!(pieces[6], [true; true;; true; true])  # bee
    push!(pieces[6], [false; true;; true; true;; true; false])
    push!(pieces[6], [false; true; false;; true; true; true])

    push!(pieces[7], [false; true; true;; false; true; false;; true; false; false])  # wave
    push!(pieces[7], [true; true; false;; false; true; false;; false; false; true])
    push!(pieces[7], [false; true; false; true;; true; false; true; false])
    push!(pieces[7], hcat([true; true; true; true]))
    push!(pieces[7], [true; false;; true; true;; false; true])
    push!(pieces[7], [false; true;; false; true;; true; false;; true; false])

    active = Dict{Pair{Int64, Int64}, Set{VariableRef}}()  #  Set of all z[w, h, p, d] that place a tile on (w, h)
    for i in 1 : width
        for j in 1 : height
            active[Pair(i, j)] = Set{VariableRef}()
        end
    end

    # Gurobi Solver
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "PoolSearchMode", 2)
    set_optimizer_attribute(model, "PoolSolutions", typemax(Int32))
    set_optimizer_attribute(model, "PoolGapAbs", 5.01)

    @variable(model, x[1:width, 1:height], Bin)  # if tile is interior
    @variable(model, y[1:width, 1:height], Bin)  # if tile is boundary
    @variable(model, z[1:width, 1:height, 1:length(pieces), 1:12], Bin)  # z[w, h, p, d] the left top corner of polyomino p in orientation d is placed on (w, h) with (1, 1) being top left conrer

    @objective(model, Max, sum(x))  # Maximize area enclosed

    """
    Constraints
    """

    # every polyomino is placed only once in one direction
    for i in 1 : length(pieces)
        @constraint(model, sum(z[:, :, i, :]) == 1)
    end

    # verify pieces dont go out of bounds and fill active set
    for w in 1 : width
        for h in 1 : height
            for p in 1 : length(pieces)
                for d in 1 : length(pieces[p])
                    piece = pieces[p][d]

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

    # if neighbour of interior is not boundary, than it is also interior: (not x[l]) or y[l'] or x[l']
    for w in 1 : width
        for h in 1 : height
            if (w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h] + x[w - 1, h] >= 1)
            end
            if (h != 1)
                @constraint(model, (1 - x[w, h]) + y[w, h - 1] + x[w, h - 1] >= 1)
            end
            if (w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h] + x[w + 1, h] >= 1)
            end
            if (h != height)
                @constraint(model, (1 - x[w, h]) + y[w, h + 1] + x[w, h + 1] >= 1)
            end

            if h % 2 == 0
                if (w != width && h != height)
                    @constraint(model, (1 - x[w, h]) + y[w + 1, h + 1] + x[w + 1, h + 1] >= 1)
                end
                if (w != width && h != 1)
                    @constraint(model, (1 - x[w, h]) + y[w + 1, h - 1] + x[w + 1, h - 1] >= 1)
                end
            else
                if (h != height && w != 1)
                    @constraint(model, (1 - x[w, h]) + y[w - 1, h + 1] + x[w - 1, h + 1] >= 1)
                end
                if (w != 1 && h != 1)
                    @constraint(model, (1 - x[w, h]) + y[w - 1, h - 1] + x[w - 1, h - 1] >= 1)
                end
            end
        end
    end

    # orientations which were left out because of symmetries are disallowed
    for p in 1 : length(pieces)
        for d in length(pieces[p]) + 1 : 12
            @constraint(model, sum(z[:, :, p, d]) == 0)
        end
    end

    # only allowed placement on parity hexagons
    for w in 1 : width
        for h in 1 : height
            if h % 2 == 0
                @constraint(model, sum(z[w, h, :, :]) == 0)
            end
        end
    end

    @constraint(model, sum(y[3, :]) >= 1)
    @constraint(model, sum(y[:, 3]) >= 1)

    optimize!(model)

    println()
    printSol(1, z, width, height, pieces)
    println()
    
    println("Number of solutions: " * string(result_count(model)))
    println()
    
    uniq = Set{Vector}()
    push!(uniq, order(1, z, width, height, pieces))

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

    writedlm("solutions-hex.csv",  output, ',')
    println("Saving completed.")
end
