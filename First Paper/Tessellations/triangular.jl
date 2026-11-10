using JuMP
using Gurobi

# https://nornagon.medium.com/til-triangle-grids-133ed47cc807

function prettyPrint(p::Matrix{Bool})
    for i in 1 : size(p, 1)
        if i % 2 == 1
            print("\\")
        else
            print("/")
        end

        for j in 1 : size(p, 2)
            if p[i, j]
                print("*")
            else
                print(" ")
            end
            if (i + j) % 2 == 1
                print("\\")
            else
                print("/")
            end
        end
        println("")
    end
end

function prettyPrint(p::Matrix{Int64})
    for i in 1 : size(p, 1)
        if i % 2 == 1
            print("\\")
        else
            print("/")
        end

        for j in 1 : size(p, 2)
            if p[i, j] != 0
                print(Char(p[i, j] - 1 + 'A'))
            else
                print(" ")
            end
            if (i + j) % 2 == 1
                print("\\")
            else
                print("/")
            end
        end
        println("")
    end
end

function solvePolyiamond(area = Pair(24, 24))
    """
    Setup
    """

    # bounding box to search in
    width = area[1]
    height = area[2]

    pieces = Vector{Vector{Matrix{Bool}}}()
    for i in 1 : 12
        push!(pieces, Vector{Matrix{Bool}}())
    end

    # hexiamond (https://mathworld.wolfram.com/Hexiamond.html)
    push!(pieces[1], [false;; true;; true;; true;; true;; true;; true])  # bar
    push!(pieces[1], [true;; true;; true;; true;; true;; true])
    push!(pieces[1], [false; false; true; true;; false; true; true; false;; true; true; false; false])
    push!(pieces[1], [false; false; true;; false; true; true;; true; true; false;; true; false; false])
    push!(pieces[1], [true; true; false; false;; false; true; true; false;; false; false; true; true])
    push!(pieces[1], [true; false; false;; true; true; false;; false; true; true;; false; false; true])

    push!(pieces[2], [false; false;; true; true;; false; true;; false; true;; false; true;; false; true])  # crook
    push!(pieces[2], [false; false;; true; true;; true; false;; true; false;; true; false;; true; false])
    push!(pieces[2], [false; false;; false; true;; false; true;; false; true;; false; true;; true; true])
    push!(pieces[2], [false; false;; true; false;; true; false;; true; false;; true; false;; true; true])
    push!(pieces[2], [false; true; true;; true; true; false;; true; false; false;; true; false; false])
    push!(pieces[2], [false; true; true;; true; true; true;; true; false; false])
    push!(pieces[2], [false; false; true;; false; false; true;; false; true; true;; true; true; false])
    push!(pieces[2], [false; false; false; false;; false; false; true; false;; true; true; true; false;; true; true; false; false])
    push!(pieces[2], [false; false; false;; true; false; false;; true; false; false;; true; true; false;; false; true; true])
    push!(pieces[2], [true; false; false;; true; true; true;; false; true; true])
    push!(pieces[2], [false; false; false;; true; true; false;; true; true; true;; false; false; true])
    push!(pieces[2], [false; false; false;; true; true; false;; false; true; true;; false; false; true;; false; false; true])

    push!(pieces[3], [false; false;; false; true;; false; true;; true; true;; false; true;; false; true])  # crown
    push!(pieces[3], [false; false;; true; false;; true; false;; true; true;; true; false;; true; false])
    push!(pieces[3], [false; true; true; false;; true; true; false; false;; true; true; false; false])
    push!(pieces[3], [false; false; false;; false; true; true;; false; true; true;; true; true; false])
    push!(pieces[3], [false; false; false;; true; true; false;; false; true; true;; false; true; true])
    push!(pieces[3], [true; true; false;; true; true; false;; false; true; true])

    push!(pieces[4], [false; true;; true; true;; false; true;; false; true;; false; true])  # sphinx
    push!(pieces[4], [false; true;; false; true;; false; true;; true; true;; false; true])
    push!(pieces[4], [true; false;; true; false;; true; false;; true; true;; true; false])
    push!(pieces[4], [true; false;; true; true;; true; false;; true; false;; true; false])
    push!(pieces[4], [false; false; false; false;; false; true; false; false;; false; true; true; false;; true; true; false; false;; true; false; false; false])
    push!(pieces[4], [true; true; true;; true; true; false;; true; false; false])
    push!(pieces[4], [false; false; false;; false; false; true;; false; true; true;; true; true; true])
    push!(pieces[4], [false; false; false;; false; false; true;; false; true; true;; true; true; true;; false; true; false])
    push!(pieces[4], [false; true; false;; true; true; false;; false; true; true;; false; false; true])
    push!(pieces[4], [false; false; false;; true; true; true;; false; true; true;; false; false; true])
    push!(pieces[4], [true; false; false;; true; true; false;; true; true; true])
    push!(pieces[4], [true; false; false;; true; true; false;; false; true; true;; false; true; false])

    push!(pieces[5], [false; true; true; false;; false; true; false; false;; false; true; false; false;; true; true; false; false])  # snake
    push!(pieces[5], [false; false; false;; true; true; false;; false; true; false;; false; true; false;; false; true; true])
    push!(pieces[5], [true; true; true;; true; true; true])
    push!(pieces[5], [false; false;; false; true;; false; true;; true; true;; true; false;; true; false])
    push!(pieces[5], [false; false;; true; false;; true; false;; true; true;; false; true;; false; true])
    push!(pieces[5], [false; false; false;; true; true; true;; true; true; true])

    push!(pieces[6], [false; true;; true; true;; false; true;; true; true])  # yacht
    push!(pieces[6], [false; false;; true; true;; false; true;; true; true;; false; true])
    push!(pieces[6], [true; false;; true; true;; true; false;; true; true])
    push!(pieces[6], [false; false;; true; true;; true; false;; true; true;; true; false])
    push!(pieces[6], [false; false; false;; false; true; false;; true; true; true;; true; true; false])
    push!(pieces[6], [false; false;; false; true;; true; true;; true; true;; true; false])
    push!(pieces[6], [false; true;; true; true;; true; true;; true; false])
    push!(pieces[6], [false; true; true;; true; true; true;; false; true; false])
    push!(pieces[6], [false; true; false; false;; true; true; true; false;; false; true; true; false])
    push!(pieces[6], [false; false;; true; false;; true; true;; true; true;; false; true])
    push!(pieces[6], [true; false;; true; true;; true; true;; false; true])
    push!(pieces[6], [false; false; false;; true; true; false;; true; true; true;; false; true; false])

    push!(pieces[7], [false; true;; false; true;; false; true;; true; true;; true; false])  # chevron
    push!(pieces[7], [true; false;; true; true;; false; true;; false; true;; false; true])
    push!(pieces[7], [false; true; true; false;; true; true; true; true])
    push!(pieces[7], [false; false; false; false;; true; true; true; true;; false; true; true; false])
    push!(pieces[7], [true; false;; true; false;; true; false;; true; true;; false; true])
    push!(pieces[7], [false; true;; true; true;; true; false;; true; false;; true; false])

    push!(pieces[8], [false; true; true;; true; true; false;; false; true; false;; false; true; false])  # signpost
    push!(pieces[8], [false; false; false;; false; true; false;; false; true; false;; true; true; false;; false; true; true])
    push!(pieces[8], [true; true;; true; true;; true; false;; true; false])
    push!(pieces[8], [true; true; true;; true; true; false;; false; true; false])
    push!(pieces[8], [false; true; false;; true; true; false;; true; true; true])
    push!(pieces[8], [false; false;; true; false;; true; false;; true; true;; true; true])
    push!(pieces[8], [false; true; false;; false; true; false;; false; true; true;; true; true; false])
    push!(pieces[8], [false; false; false;; true; true; false;; false; true; true;; false; true; false;; false; true; false])
    push!(pieces[8], [true; true;; true; true;; true; false;; true; false])
    push!(pieces[8], [false; false; false;; false; true; false;; false; true; true;; true; true; true])
    push!(pieces[8], [false; false; false;; true; true; true;; false; true; true;; false; true; false])
    push!(pieces[8], [true; true;; true; true;; false; true;; false; true])

    push!(pieces[9], [false; true; true;; true; true; false;; false; true; true])  # lobster
    push!(pieces[9], [false; false; false;; true; true; false;; false; true; true;; true; true; false])
    push!(pieces[9], [false; false;; false; true;; true; true;; true; true;; false; true])
    push!(pieces[9], [false; true;; true; true;; true; true;; false; true])
    push!(pieces[9], [false; false;; true; false;; true; true;; true; true;; true; false])
    push!(pieces[9], [true; false;; true; true;; true; true;; true; false])

    push!(pieces[10], [false; true;; true; true;; true; false;; true; true])  # hook
    push!(pieces[10], [false; false;; true; true;; true; false;; true; true;; false; true])
    push!(pieces[10], [true; false;; true; true;; false; true;; true; true])
    push!(pieces[10], [false; false;; true; true;; false; true;; true; true;; true; false])
    push!(pieces[10], [false; true;; false; true;; true; true;; true; true])
    push!(pieces[10], [false; false; true;; true; true; true;; false; true; true])
    push!(pieces[10], [false; false; false;; true; true; false;; true; true; true;; true; false; false])
    push!(pieces[10], [false; false;; true; true;; true; true;; true; false;; true; false])
    push!(pieces[10], [false; true; true;; true; true; true;; false; false; true])
    push!(pieces[10], [false; false;; true; true;; true; true;; false; true;; false; true])
    push!(pieces[10], [true; false;; true; false;; true; true;; true; true])
    push!(pieces[10], [false; false; false;; true; false; false;; true; true; true;; true; true; false])

    push!(pieces[11], [false; false;; true; true;; true; true;; true; true])  # hexagon

    push!(pieces[12], [true; true;; true; true;; true; true])  # butterfly
    push!(pieces[12], [false; false; false;; false; true; false;; false; true; true;; true; true; false;; false; true; false])
    push!(pieces[12], [false; true; false;; true; true; false;; false; true; true;; false; true; false])

    active = Dict{Pair{Int64, Int64}, Set{VariableRef}}()  #  Set of all z[w, h, p, d] that place a tile on (w, h)
    for i in 1 : width
        for j in 1 : height
            active[Pair(i, j)] = Set{VariableRef}()
        end
    end

    # Gurobi Solver
    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:width, 1:height], Bin)  # if tile is interior
    @variable(model, y[1:width, 1:height], Bin)  # if tile is boundary
    @variable(model, z[1:width, 1:height, 1:length(pieces), 1:12], Bin)  # z[w, h, p, d] the left top corner of polyomino p in orientation d is placed on (w, h) with (1, 1) being top left conrer

    @objective(model, Max, sum(x))  # Maximize area enclosed

    """
    Constraints
    """

    # every polyomino is placed only once in one direction
    for i in 1 : 12
        @constraint(model, sum(z[:, :, i, :]) == 1)
    end

    # verify pieces dont go out of bounds and fill active set
    for w in 1 : width
        for h in 1 : height
            for p in 1 : 12
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
            if (w + h) % 2 == 1
                if (h != 1)
                    @constraint(model, (1 - x[w, h]) + y[w, h - 1] + x[w, h - 1] >= 1)  # above
                end
                if (w != width)
                    @constraint(model, (1 - x[w, h]) + y[w + 1, h] + x[w + 1, h] >= 1)  # right
                end
                if (h != height && w != width && w != width - 1)
                    @constraint(model, (1 - x[w, h]) + y[w + 2, h + 1] + x[w + 2, h + 1] >= 1)  # right below
                end
                if (h != height)
                    @constraint(model, (1 - x[w, h]) + y[w, h + 1] + x[w, h + 1] >= 1)  # below
                end
                if (h != height && w != 1 && w != 2)
                    @constraint(model, (1 - x[w, h]) + y[w - 2, h + 1] + x[w - 2, h + 1] >= 1)  # left below
                end
                if (w != 1)
                    @constraint(model, (1 - x[w, h]) + y[w - 1, h] + x[w - 1, h] >= 1)  # left
                end
            else
                if (h != 1)
                    @constraint(model, (1 - x[w, h]) + y[w, h - 1] + x[w, h - 1] >= 1)  # above
                end
                if (h != 1 && w != width && w != width - 1)
                    @constraint(model, (1 - x[w, h]) + y[w + 2, h - 1] + x[w + 2, h - 1] >= 1)  # right above
                end
                if (w != width)
                    @constraint(model, (1 - x[w, h]) + y[w + 1, h] + x[w + 1, h] >= 1)  # right
                end
                if (h != height)
                    @constraint(model, (1 - x[w, h]) + y[w, h + 1] + x[w, h + 1] >= 1)  # below
                end
                if (w != 1)
                    @constraint(model, (1 - x[w, h]) + y[w - 1, h] + x[w - 1, h] >= 1)  # left
                end
                if (h != 1 && w != 1 && w != 2)
                    @constraint(model, (1 - x[w, h]) + y[w - 2, h - 1] + x[w - 2, h - 1] >= 1)  # left above
                end
            end

            # middle neighbours (connected 120 degrees)
            if (w != 1 && h != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h - 1] + x[w - 1, h - 1] >= 1)  # left above
            end
            if (h != 1 && w != width)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h - 1] + x[w + 1, h - 1] >= 1)  # right above
            end
            if (h != height && w != 1)
                @constraint(model, (1 - x[w, h]) + y[w - 1, h + 1] + x[w - 1, h + 1] >= 1)  # left below
            end
            if (w != width && h != height)
                @constraint(model, (1 - x[w, h]) + y[w + 1, h + 1] + x[w + 1, h + 1] >= 1)  # right below
            end
            if (w != 1 && w != 2)
                @constraint(model, (1 - x[w, h]) + y[w - 2, h] + x[w - 2, h] >= 1)  # left above
            end
            if (w != width && w != width - 1)
                @constraint(model, (1 - x[w, h]) + y[w + 2, h] + x[w + 2, h] >= 1)  # right below
            end
        end
    end

    for p in 1 : 12  # orientations which were left out because of symmetries are disallowed
        for d in length(pieces[p]) + 1 : 12
            @constraint(model, sum(z[:, :, p, d]) == 0)
        end
    end

    # only allowed placement on downwards pointing triangles
    for w in 1 : width
        for h in 1 : height
            if (w + h) % 2 == 1
                @constraint(model, sum(z[w, h, :, :]) == 0)
            end
        end
    end

    optimize!(model)

    output = fill(0, width, height)

    solution = findall(t->t == 1.0, round.(value.(z)))
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

    display(output)

    println()

    prettyPrint(output)
end
