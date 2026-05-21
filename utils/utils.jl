module Utils

    export refineGrid1d

    function refineGrid1d(grid::Vector{Float64})
        grid_mid = (grid[1:end-1] .+ grid[2:end])/2
        return sort(vcat(grid, grid_mid))
    end

end