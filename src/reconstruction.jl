function initial_guess(size, p :: AbstractFloat)
    array   = zeros(Bool, size)
    nsolid  = (1 - p)*length(array) |> round |> Int
    n       = 0
    indices = CartesianIndices(array)

    while true
        idx = rand(indices)

        if array[idx] == false
            n = n + 1
            array[idx] = true
        end

        if n == nsolid
            return array
        end
    end
end

function porosity(s2ft :: AbstractArray{<: AbstractFloat}, size)
    s2 = irfft(s2ft, size[1]) / reduce(*, size)
    return 1 - s2[1]
end

"""
    two_point(array)

Calculate unnormalized two-point correlation function for array of
booleans. The result can be used in [`phaserec`](@ref) function.
"""
two_point(array :: AbstractArray{Bool}) = array |> rfft .|> abs2

# This is the core of PhaseRec algorithm.
#
# This function takes an array and replaces absolute value of its FFT
# with content of s2ft and then applies low-pass filter to remove
# noise.
function replace_abs(array   :: AbstractArray{Bool, N},
                     lowpass :: AbstractArray{C, N},
                     s2ft    :: AbstractArray{R, N}) where {C <: Complex, R <: AbstractFloat, N}
    myifft(x) = irfft(x, size(array, 1))
    ft = rfft(array)

    repft = @. ft * sqrt(s2ft) / abs(ft)
    repft = replace(x -> isnan(x) ? 0 : x, repft)
    repft = repft .* lowpass

    return repft |> myifft .|> real
end

# Make Gaussian low-pass filter
function make_filter(size, σ)
    array = zeros(Float64, size)
    w = 4ceil(Int, σ) + 1
    for idx in CartesianIndices(array)
        x = ((x-1-w)^2 for x in Tuple(idx))
        array[idx] = exp(-sum(x)/(2σ^2))
    end

    array = array ./ sum(array)
    return rfft(array)
end

# Convert grayscale array to binary saving porosity of the original
threshold(array :: AbstractArray{<: AbstractFloat}, porosity :: AbstractFloat) =
    array .> quantile(reshape(array, length(array)), porosity)

"""
    phaserec(s2ft, size; radius = 0.6, maxstep = 300, ϵ = 10^-5)

Reconstruct an image from the two-point correlation function. `s2ft`
is unnormalized two-point correlation function in frequency domain and
`size` are dimensions of the original image. Optional parameter
`radius` governs noise filtration. `maxsteps` is a maximal number of
iterations. Smaller `ϵ` usually produces better results, but default
value should be OK for all cases.

See also: [`two_point`](@ref).
"""
function phaserec(s2ft :: AbstractArray{<: AbstractFloat}, size;
                  radius = 0.6, maxsteps = 300, ϵ = 10^-5)
    # Make initial guess of
    p = porosity(s2ft, size)
    recon  = initial_guess(size, p)
    filter = make_filter(size, radius)
    normfn(corr, rec) = norm((corr - two_point(rec))/length(corr))
    initnorm = normfn(s2ft, recon)
    oldn = 1.0

    for steps in 1:maxsteps
        gray = replace_abs(recon, filter, s2ft)
        recon = threshold(gray, p)

        n = normfn(s2ft, recon) / initnorm

        if (rem(steps, 10) == 1)
            @info "Cost = $(n)"
        end

        if oldn - n < ϵ || n > oldn
            break
        end
        oldn = n
    end

    return recon
end

"""
    phaserec(array; radius = 0.6, maxsteps = 300, ϵ = 10^-5)

Reconstruct `array` which must be an array of booleans.
"""
phaserec(array :: AbstractArray{Bool}; radius = 0.6, maxsteps = 300, ϵ = 10^-5) =
    phaserec(two_point(array), size(array);
             radius = radius, maxsteps = maxsteps, ϵ = ϵ)
