struct Pentomino
    id::Int
    name::Char
    shape::Vector{Tuple{Int, Int}}  # list of (x, y) coordinates
end

struct Possibility
    new_orient::Int
    area_change_8x::Int
    pos_diff::Tuple{Int, Int}

    # data for reconstructing the pentomino placement
    pento_name::Char
    shifted_tiles::Vector{Tuple{Int, Int}}
    next_tile_offset::Tuple{Int, Int}
end

const raw_pentos = [
    Pentomino(1, 'F', [(0, 1), (1, 1), (1, 0), (1, 2), (2, 0)]),
    Pentomino(2, 'I', [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4)]),
    Pentomino(3, 'L', [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0)]),
    Pentomino(4, 'N', [(0, 0), (0, 1), (1, 1), (1, 2), (1, 3)]),
    Pentomino(5, 'P', [(0, 0), (0, 1), (1, 0), (1, 1), (0, 2)]),
    Pentomino(6, 'T', [(0, 0), (1, 0), (2, 0), (1, 1), (1, 2)]),
    Pentomino(7, 'U', [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1)]),
    Pentomino(8, 'V', [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]),
    Pentomino(9, 'W', [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)]),
    Pentomino(10, 'X', [(1, 0), (0, 1), (1, 1), (2, 1), (1, 2)]),
    Pentomino(11, 'Y', [(0, 0), (0, 1), (0, 2), (0, 3), (1, 1)]),
    Pentomino(12, 'Z', [(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)])
]

function normalize(P)
    vals_x = [x for (x, y) in P]
    vals_y = [y for (x, y) in P]
    min_x = minimum(vals_x)
    min_y = minimum(vals_y)
    
    return sort([(x - min_x, y - min_y) for (x, y) in P])  # sort lexicographically
end

function rotate(P, num_rot)
    if num_rot == 0
        return P
    elseif num_rot == 1
        return [(-y, x) for (x, y) in P]
    elseif num_rot == 2
        return [(-x, -y) for (x, y) in P]
    else
        return [(y, -x) for (x, y) in P]
    end
end

function flip(P)
    return [(y, x) for (x, y) in P]
end

function get_orientations(P)
    base = [normalize(rotate(P, r)) for r in 0:3]
    flipped = [normalize(rotate(flip(P), r)) for r in 0:3]
    return unique(vcat(base, flipped))
end

# Normalize all given pentominoes
const pentos = [Pentomino(p.id, p.name, normalize(p.shape)) for p in raw_pentos]

function pretty_print(P)
    P = normalize(P) # Normalize before drawing
    xs = first.(P)
    ys = last.(P)

    xmax = maximum(xs)
    ymax = maximum(ys)

    grid = fill('.', ymax + 1, xmax + 1)  # range from 0 to max, so we need max + 1 much space in each dimension

    for (x,y) in P
        grid[ymax - y + 1, x + 1] = '#'  # print in cartesian coordinates
    end

    # Print row by row
    for row in eachrow(grid)
        println(join(row))
    end
end

for p in pentos
    orients = get_orientations(p.shape)

    for (i, shape) in enumerate(orients)
        println("Pentomino $(p.name), Orientation $i:")
        pretty_print(shape)
        println()  # extra visual spacing before next orientation
    end
end

const DIRECTIONS = [(0, -1), (1, 0), (0, 1), (-1, 0)]

function compute_boundary(P)
    boundary_successor = Dict{Tuple{Int, Int}, Tuple{Int, Int}}()  # key: coordinates of starting point, value: coordinates of ending point
    edge_to_tile = Dict{Tuple{Tuple{Int, Int}, Tuple{Int, Int}}, Tuple{Tuple{Int, Int}, Tuple{Int, Int}}}()  # (perimeter edge) -> (adjacent inner tile, adjacent outer tile)

    for (x, y) in P
        for (dx, dy) in DIRECTIONS
            xx, yy = x + dx, y + dy
            if !((xx, yy) in P)  # pentomino tile (x, y) shares outside edge with non-pentomino tile (xx, yy)
                
                # (dx, dy) is a vector pointing from (x, y) to (xx, yy)
                # then (-dy, dx) is a vector rotated 90 deg counterclockwise (since (-dy, dx) = M * (dx, dy) where M is the 90 deg ccw rotation from here: https://en.wikipedia.org/wiki/Rotation_matrix#Common_2D_rotations)
                perpx2, perpy2 = -dy, dx  # half length to allow easier start and midpoint creation
                
                midx2, midy2 = x + xx, y + yy  # midpoint of the edge
                
                # calculate start and end points of red attachement vectors
                start_pt = (midx2 + perpx2, midy2 + perpy2)
                end_pt = (midx2 - perpx2, midy2 - perpy2)
                
                # store the clockwise successor relationship
                @assert !haskey(boundary_successor, start_pt)
                boundary_successor[start_pt] = end_pt

                edge_to_tile[(start_pt, end_pt)] = ((x, y), (xx, yy))
            end
        end
    end

    return boundary_successor, edge_to_tile
end

function construct_boundary_path(P)
    boundary_successor, edge_to_tile = compute_boundary(P)
    
    # construct clockwise edge path as discussed
    start_pt = current = minimum(keys(boundary_successor))  # can be any edge, but for debugging choose deterministically
    path = [start_pt]
    
    current = boundary_successor[current]
    while current != start_pt
        push!(path, current)
        current = boundary_successor[current]
    end
    
    @assert length(path) == length(boundary_successor)
    
    return path, edge_to_tile
end

const DIR_TO_INDEX = Dict(
    (0, -2) => 1,
    (2,  0) => 2,
    (0,  2) => 3,
    (-2, 0) => 4,
)

const possibilities = Dict{Tuple{Int, Int}, Vector{Possibility}}()  # (pento, start_orient) -> (new_orient, area_change, pento_size) 

for p in pentos
    orients = get_orientations(p.shape)  # get all orientations (including flipping)
    # unique_orients handled by get_orientations returning unique

    for (j, tiles) in enumerate(orients)
        path, edge_to_tile = construct_boundary_path(tiles)
        doubled_path = vcat(path, path)  # allow wrap-around

        len_path = length(path)
        for start_index in 1:len_path
            for end_index in (start_index + 1):(start_index + len_path - 1)
                
                # denote orientation as difference vector between two endpoints of attachment edge
                p1 = doubled_path[start_index]
                p2 = doubled_path[start_index + 1]
                diff_start = (p2[1] - p1[1], p2[2] - p1[2])
                start_orient = DIR_TO_INDEX[diff_start]

                # reversed since new S is E + 1 and new S + 1 is E
                p3 = doubled_path[end_index]
                p4 = doubled_path[end_index + 1]
                diff_end = (p3[1] - p4[1], p3[2] - p4[2])
                new_orient = DIR_TO_INDEX[diff_end]

                sols = get!(possibilities, (p.id, start_orient), Possibility[])

                # extract path of inner edges translated to S + 1 at (0, 0)
                origin = doubled_path[start_index + 1]
                inner_path = [
                    (pt[1] - origin[1], pt[2] - origin[2])
                    for pt in doubled_path[(start_index + 1):end_index]
                ]

                # pre-compute area_change with shoelace formula
                area_change_8x = 0
                for e in 1:(length(inner_path) - 1)
                    area_change_8x += (
                        inner_path[e][1] * inner_path[e + 1][2]
                        - inner_path[e + 1][1] * inner_path[e][2]
                    )  # shoelace term (8x since 1/2 is missing and both coordinates are multiplied by 2x)
                end

                # calculate data for reconstruction

                start_tile_key = (doubled_path[start_index], doubled_path[start_index + 1])
                start_tile, _ = edge_to_tile[start_tile_key]  # inner tile of edge
                shifted_tiles = [(t[1] - start_tile[1], t[2] - start_tile[2]) for t in tiles]  # shift tiles such that start_tile has coordinates (0, 0)
                
                end_tile_key = (doubled_path[end_index], doubled_path[end_index + 1])
                _, next_tile = edge_to_tile[end_tile_key]  # outer tile of edge
                shifted_next_tile = (next_tile[1] - start_tile[1], next_tile[2] - start_tile[2])

                poss = Possibility(
                    new_orient,
                    area_change_8x,
                    inner_path[end],
                    p.name,
                    shifted_tiles,
                    shifted_next_tile
                )

                push!(sols, poss)  # inner_path[0] is always (0, 0), so don't need to save
            end
        end
    end
end

max_poss = maximum(length(possib) for ((_, _), possib) in possibilities)
println("No combination has more than $max_poss possibilities.")

# State packing
@inline function pack_state(pos, orient)
    # x: 12 bits, y: 12 bits, orient: 4 bits

    # shift from [-2048, 2047] to [0, 4095] to prevent negative values
    p_x = Int32(pos[1] + 2048)
    p_y = Int32(pos[2] + 2048)
    return (p_x << 16) | (p_y << 4) | Int32(orient)
end

@inline function unpack_state(packed)
    orient = Int(packed & 15)
    y = Int((packed >> 4) & 4095) - 2048
    x = Int((packed >> 16) & 4095) - 2048
    return ((x, y), orient)
end

@inline function count_set_bits(n::UInt16)
    c = 0
    while n > 0
        n &= (n - 1)
        c += 1
    end
    return c
end

@inline function update_solutions!(new_solutions, new_remember, state_key, new_area_8x, info, new_pos, new_orient, start_pos, start_orient, size, max_fence_area, optimal_subsets, new_used, no_output=false)
    found = haskey(new_solutions, state_key)

    if !found || new_area_8x > new_solutions[state_key]
        new_solutions[state_key] = new_area_8x
        new_remember[state_key] = info

        if new_pos == start_pos && new_orient == start_orient
            @assert new_area_8x % 8 == 0  # area must be integer

            true_area = div(new_area_8x, 8)
            if true_area > max_fence_area
                max_fence_area = true_area
                if !no_output
                    println("Size: $size, Max Fence Area: $max_fence_area")
                    flush(stdout)  # ensure immediate output
                end

                optimal_subsets = new_used
            end
        end
    end
    return max_fence_area, optimal_subsets  # since optimal_subsets is now just a number, we have to return it explicitly
end

@inline function update_solutions_all!(new_solutions, new_remember, state_key, new_area_8x, info, new_pos, new_orient, start_pos, start_orient, size, max_fence_area, optimal_subsets, new_used, no_output=false)
    found = haskey(new_solutions, state_key)

    if !found || new_area_8x >= new_solutions[state_key]  # either the first time we reach this position with those pentominoes or we improve the area
        if found && new_area_8x == new_solutions[state_key]
            push!(new_remember[state_key], info)  # equal to current best solution, collect as another predecessors
        else
            new_solutions[state_key] = new_area_8x
            new_remember[state_key] = [info]  # better than all previous values, reset collection of predecessors
        end

        if new_pos == start_pos && new_orient == start_orient
            @assert new_area_8x % 8 == 0  # area must be integer

            true_area = div(new_area_8x, 8)
            if true_area > max_fence_area
                max_fence_area = true_area
                if !no_output
                    println("Size: $size, Max Fence Area: $max_fence_area")
                    flush(stdout)  # ensure immediate output
                end

                empty!(optimal_subsets)
                push!(optimal_subsets, new_used)
            elseif true_area == max_fence_area
                push!(optimal_subsets, new_used)
            end
        end
    end
    return max_fence_area
end

struct ReconstructionStep
    name::Char
    tiles::Vector{Tuple{Int, Int}}
    next_tile::Tuple{Int, Int}
end

function find_single_path(used_mask, remember)
    packed_state = pack_state((0, 0), 3)
    used = used_mask

    path = Vector{ReconstructionStep}()

    while used != 0
        info = remember[used][packed_state]  # find which pentomino was used to reach current state

        # Unpack info: idx << 7 | start_orient << 4 | pento_idx
        idx = (info >> 7)
        orient = (info >> 4) & 0x07
        pento_idx = info & 0x0F

        poss = possibilities[(pento_idx, orient)][idx]
        
        # Reconstruct previous state
        (curr_pos, current_end_orient) = unpack_state(packed_state)
        prev_pos = (curr_pos[1] - poss.pos_diff[1], curr_pos[2] - poss.pos_diff[2])
        packed_state = pack_state(prev_pos, orient)

        used = used & ~(UInt16(1) << (pento_idx - 1))

        push!(path, ReconstructionStep(poss.pento_name, poss.shifted_tiles, poss.next_tile_offset))
    end

    return reverse(path)
end

function find_all_paths!(all_paths, used_mask, remember)
    # Stack: (packed_state, used, path)
    stack = Any[(pack_state((0, 0), 3), used_mask, ReconstructionStep[])]
    
    while !isempty(stack)
        (packed_state, used, path) = pop!(stack)

        if used == 0
            push!(all_paths, reverse(path))
            continue
        end
        
        for info in remember[used][packed_state]
            # Unpack info: idx << 7 | start_orient << 4 | pento_idx
            idx = (info >> 7)
            orient = (info >> 4) & 0x07
            pento_idx = info & 0x0F

            # Get possibility again
            poss = possibilities[(pento_idx, orient)][idx]
            
            # Reconstruct previous state
            (curr_pos, current_end_orient) = unpack_state(packed_state)
            # Reverse the position change
            prev_pos = (curr_pos[1] - poss.pos_diff[1], curr_pos[2] - poss.pos_diff[2])
            # The previous orientation was the start_orient of this move
            prev_packed = pack_state(prev_pos, orient)

            prev_used = used & ~(UInt16(1) << (pento_idx - 1))
            
            new_path = copy(path)
            push!(new_path, ReconstructionStep(poss.pento_name, poss.shifted_tiles, poss.next_tile_offset))
            
            push!(stack, (prev_packed, prev_used, new_path))
        end
    end
end

function process_path(path)
    fence_tiles = Dict{Tuple{Int, Int}, Char}()
    curr = (0, 0)
    
    for step in path
        for (x, y) in step.tiles
            @assert !haskey(fence_tiles, (x + curr[1], y + curr[2]))  # Check for self-Intersecting pentominoes
            fence_tiles[(x + curr[1], y + curr[2])] = step.name
        end
        curr = (curr[1] + step.next_tile[1], curr[2] + step.next_tile[2])  # update current coordinates 
    end
    return fence_tiles
end

# Print the grid visualization
function print_grid(fence_tiles; io::IO = stdout)
    if isempty(fence_tiles)
        println("Empty grid")
        return
    end
    
    xs = [k[1] for k in keys(fence_tiles)]
    ys = [k[2] for k in keys(fence_tiles)]
    min_x, max_x = minimum(xs), maximum(xs)
    min_y, max_y = minimum(ys), maximum(ys)
    
    grid = [['.' for _ in min_x:max_x] for _ in min_y:max_y]
    
    for ((x, y), c) in fence_tiles
        grid[max_y - y + 1][x - min_x + 1] = c
    end
    
    for row in grid
        println(io, join(row))
    end
end

function calculate_enclosed_areas(fence_tiles)
    if isempty(fence_tiles)
        return Int[]
    end

    xs = [k[1] for k in keys(fence_tiles)]
    ys = [k[2] for k in keys(fence_tiles)]
    min_x, max_x = minimum(xs), maximum(xs)
    min_y, max_y = minimum(ys), maximum(ys)
    
    areas = Int[]
    visited = Set{Tuple{Int, Int}}()  # Track visited cells
    dirs8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]  # 8-connectivity: all 8 directions including diagonals
    
    for x in min_x:max_x
        for y in min_y:max_y
            if !haskey(fence_tiles, (x, y)) && !((x, y) in visited)
                # BFS search for all connected components in the complement that are not connected to the outside
                push!(visited, (x, y))
                q = [(x, y)]
                size = 1
                touches_border = false
                
                head = 1
                while head <= length(q)
                    cx, cy = q[head]
                    head += 1
                    
                    if cx == min_x || cx == max_x || cy == min_y || cy == max_y
                        touches_border = true
                    end
                    
                    for (dx, dy) in dirs8
                        nx, ny = cx + dx, cy + dy
                        if nx >= min_x && nx <= max_x && ny >= min_y && ny <= max_y
                            if !haskey(fence_tiles, (nx, ny)) && !((nx, ny) in visited)
                                push!(visited, (nx, ny))
                                push!(q, (nx, ny))
                                size += 1
                            end
                        end
                    end
                end
                if !touches_border
                    push!(areas, size)
                end
            end
        end
    end
    sort!(areas, rev=true)  # sort in descending order
    return areas
end

function optimal_order_fence!(order, table = nothing, remember = nothing; find_all_solutions = false)
    num_pieces = length(order)

    # Precompute masks for all prefixes
    masks = Vector{UInt16}(undef, num_pieces + 1)
    masks[1] = UInt16(0)
    current_mask = UInt16(0)
    for i in 1:num_pieces
        current_mask |= (UInt16(1) << (order[i] - 1))
        masks[i+1] = current_mask
    end

    if remember === nothing
        # use UInt16 to store (idx << 7) | (start_orient << 4) | pento_idx
        remember = find_all_solutions ? Dict{UInt16, Dict{Int32, Vector{UInt16}}}() : Dict{UInt16, Dict{Int32, UInt16}}()
    end

    optimal_subsets = find_all_solutions ? Set{UInt16}() : UInt16(0)  # not strictly necessary, but kept for compatibility 

    max_fence_area = 0

    start_pos = (0, 0)
    start_orient = 3
    start = pack_state(start_pos, start_orient)

    start_index = 0

    if table === nothing || length(table) != num_pieces + 1
        table = Vector{Dict{UInt16, Dict{Int32, Int}}}(undef, num_pieces + 1)  # map set of placed pentominoes to all possible placements
    else
        # findest longest prefix of 'order' that is already computed
        while start_index < num_pieces - 1 && haskey(table[start_index + 2], masks[start_index + 2])
            start_index += 1
        end
    end

    if start_index == 0
        table[1] = Dict(UInt16(0) => Dict(start => 0))  # for no pentomino placed (empty set) there is only the solution that we are still at (0, 0) looking east and have 0 partial area
    end

    for size in (start_index + 1):num_pieces
        table[size+1] = Dict{UInt16, Dict{Int32, Int}}()
        
        pento_idx = order[size]

        for (used, solutions) in table[size]
            new_used = masks[size+1]  # simply choose order prefix

            # collect solutions for new used set, new_solutions is just a reference to table[size+1][new_used]
            new_solutions = get!(table[size+1], new_used) do  # do is performance improvement, since dict is only built when necessary
                Dict{Int32, Int}()
            end 

            new_remember = get!(remember, new_used) do
                find_all_solutions ? Dict{Int32, Vector{UInt16}}() : Dict{Int32, UInt16}()
            end

            for (packed_state, area_8x) in solutions
                (pos, orient) = unpack_state(packed_state)
                
                for (idx, poss) in enumerate(possibilities[(pento_idx, orient)])
                    trans_corr = pos[2] * (-poss.pos_diff[1]) + pos[1] * (poss.pos_diff[2])  # correction for precomputed value since we have a translation from (0, 0) to pos
                    new_area_8x = area_8x + poss.area_change_8x + trans_corr

                    new_pos = pos .+ poss.pos_diff
                    
                    if size == num_pieces && !(new_pos == start_pos && poss.new_orient == start_orient)
                        continue  # all pentominoes are placed but fence is not closed
                    end

                    info = UInt16((idx << 7) | (orient << 4) | pento_idx)  # data for reconstruction
                    
                    state_key = pack_state(new_pos, poss.new_orient)

                    if !find_all_solutions
                        max_fence_area, optimal_subsets = update_solutions!(new_solutions, new_remember, state_key, new_area_8x, info, new_pos, poss.new_orient, start_pos, start_orient, size, max_fence_area, optimal_subsets, new_used, true)
                    else
                        max_fence_area = update_solutions_all!(new_solutions, new_remember, state_key, new_area_8x, info, new_pos, poss.new_orient, start_pos, start_orient, size, max_fence_area, optimal_subsets, new_used, true)
                    end
                end
            end
        end
    end

    if max_fence_area == 0
        println("No solution found.")
        return
    end

    if !find_all_solutions
        sol_path = find_single_path(optimal_subsets, remember)
        return max_fence_area, sol_path, table, remember
    end

    all_paths = Vector{Vector{ReconstructionStep}}()  # collect all solution paths

    for used in optimal_subsets
        find_all_paths!(all_paths, used, remember)  # find all paths for each optimal subset of pentominoes
    end

    return max_fence_area, all_paths, table, remember
end

using Combinatorics
using ProgressMeter

function search()
    if length(ARGS) < 2
        println("Usage: julia pentomino-minmax.jl <first> <third>")
        exit(1)
    end

    first = parse(Int, ARGS[1])
    third = parse(Int, ARGS[2])
    second = 1

    mkpath("data")
    output_file = "data/order-$first-$third.csv"

    table = nothing
    remember = nothing
    
    # Progress bar
    p = Progress(factorial(9), desc="Searching permutations: ")
    
    # Open file for writing
    open(output_file, "w") do f
        remain = [n for n in 1 : 12 if !(n in (first, second, third))]

        for partial in permutations(remain)
            order = [first, second, third, partial...]

            max_fence_area, path, table, remember = optimal_order_fence!(order, table, remember)
            
            fence_tiles = process_path(path)  # verifies that no self-intersection occurs
            areas = calculate_enclosed_areas(fence_tiles)
            @assert length(areas) == 1 && areas[1] == max_fence_area
            
            println(f, "$max_fence_area,$(join(order, ','))")
            flush(f)
            
            next!(p)  # update progress bar
        end
    end
    
    finish!(p)
end

search()