module SpecUtils

using StatsBase, PyPlot, Statistics

export combine, binning


"""
    `combine(x::Array{Array{T,1},1}, y::Array{Array{T,N},1}) where {T, N}`

Combines and sorts datasets. x values of each set is stored in `x` y values
are stored in `y`.
"""
function combine(x::Array{Array{T,1},1}, y::Array{Array{T,N},1}) where {T, N}
    xc = vcat(x...)
    yc = vcat(y...)
    xyc = hcat(xc, yc)
    xys = sortslices(xyc, dims=1)
    xc .= xys[:,1]
    N == 2 && (yc .= xys[:,2:end])
    N == 1 && (yc .= xys[:,2])
    return xc, yc
end


"""
    `binning(x::Array{T,1}, y::Array{T,1}, n::Int; mode=:in, debug=false) where T <: Number`

Averages the array `y` appropriately over a certain amount of 'x' values
to return an array with 'n' data points.

Keyword Arguments:
    * `binning` = `:in` or `:out` controls the position of the first
        and last bin. `:out` centers the first and last data point in the first and
        last bin respectively. `:in` places the first bins left edge to the first
        data point and the last bins right edge to the last data point.

    * `debug` = `true` or `false`: Generate a plot for debugging.
"""
function binning(x::Array{T,1}, y::Array{T,1}, n::Int; mode=:in, debug=false) where T <: Number
    @assert length(x) == length(y)
    @assert issorted(x)

    xr = Array{Float64,1}(undef, n)
    yr = Array{Float64,1}(undef, n)

    # We want to place the edges in such a way that the data points lie exactly
    # in the middle of the bins when the number of resampled points matches the
    # number of native points. Therefore edges_tight has to be adjusted.
    weight = ProbabilityWeights(ones(T, length(x)))

    if mode == :out
        edges = range(x[1], stop=x[end], length=n)
        edges_start = x[1]   - 1/2 * step(edges)
        edges_end   = x[end] + 1/2 * step(edges)
        edges = range(edges_start, edges_end, length=n+1)
    elseif mode == :in
        edges = range(x[1] - eps(Float64), stop=x[end], length=n+1)
    else
        error("Mode $mode not defined. Choose :in or :out.")
    end

    hist = fit(Histogram, x, weight, edges, closed=:right)

    i = 1
    j = 1
    while j <= n
        w = trunc(Int, i + hist.weights[j] - 1)
        yr[j] = mean(y[i:w])
        xr[j] = (edges[j] + edges[j+1]) / 2
        i += trunc(Int, hist.weights[j])
        j += 1
    end

    if debug == true
        figure()
        title("Binning Debug")
        scatter(x, y)
        plot(xr, yr, "r.-")
        vlines(edges, ylim()..., alpha=0.5, color="k")
        for i = 1:length(hist.weights)
            isodd(i) ? (p = 0.85) : (p = 0.81)
            text(xr[i], p * ylim()[2] - ylim()[1], trunc.(Int, hist.weights[i]),
                horizontalalignment="center",
                color="k")
        end
        show()
    end

    xr, yr
end

end # module
