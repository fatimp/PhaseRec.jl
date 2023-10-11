Maybe{T} = Union{T, Nothing}

abstract type BoundaryConditions end
struct ZeroPadding <: BoundaryConditions end
struct Periodic    <: BoundaryConditions end

maybe_pad(array :: AbstractArray{Bool}, :: Periodic) = array
function maybe_pad(array :: AbstractArray{Bool, N}, :: ZeroPadding) where N
    newsize = Tuple(2s - 1 for s in size(array)) :: NTuple{N, Int64}
    indices = Tuple(Base.OneTo(s) for s in size(array)) :: NTuple{N, Base.OneTo}
    padded = falses(newsize...)
    padded[indices...] .= array
    return padded
end

maybe_cut_padding(array :: AbstractArray, :: Periodic) = array
function maybe_cut_padding(array :: AbstractArray, :: ZeroPadding)
    newsize = (Base.OneTo((s + 1)÷2) for s in size(array))
    return array[newsize...]
end

# Convert grayscale array to binary saving porosity of the original
threshold(array :: AbstractArray{<: AbstractFloat}, porosity :: AbstractFloat) =
    array .> quantile(reshape(array, length(array)), porosity)

initial_guess(size, p :: AbstractFloat) =
    threshold(rand(Float64, size), p)

function porosity(s2ft :: AbstractArray{<: AbstractFloat}, size, bc)
    recsize = bc == Periodic() ? size[1] : 2size[1] - 1
    s2 = irfft(s2ft, recsize)
    return 1 - s2[1] / reduce(*, size)
end

"""
    autocorrelation(array)

Calculate unnormalized autocorrelationcorrelation for array of
booleans. The result can be used in [`phaserec`](@ref) function.
"""
autocorrelation(array :: AbstractArray{Bool}) = array |> rfft .|> abs2

# This is the core of PhaseRec algorithm.
#
# This function takes an array and replaces absolute value of its FFT
# with content of s2ft and then does inverse FFT.
function replace_abs(array   :: AbstractArray{Bool, N},
                     s2ft    :: AbstractArray{R, N}) where {R <: AbstractFloat, N}
    ft = rfft(array)
    repft = @. ft * sqrt(s2ft) / abs(ft)

    return replace(x -> isnan(x) ? 0 : x, repft)
end

# Make Gaussian low-pass filter
function make_filter(size, σ)
    array = zeros(Float64, size)
    # Filter width
    l = 4ceil(Int, σ) + 1
    w = l >> 1

    uidx = array |> CartesianIndices |> first |> oneunit
    for offset in -w*uidx:w*uidx
        k = Tuple(offset)
        idx = (mod(k + s, s) + 1 for (k, s) in zip(k, size)) |>
            Tuple |> CartesianIndex

        array[idx] = exp(-sum(k^2 for k in k)/(2σ^2))
    end

    return rfft(array ./ sum(array))
end

"""
    phaserec(s2ft, size; radius = 0.6, maxstep = 300, ϵ = 10^-5[, noise])

Reconstruct an image from its autocorrelation. `s2ft` is unnormalized
autocorrelation in frequency domain and `size` are dimensions of the
original image. Optional parameter `radius` governs noise
filtration. `maxsteps` is a maximal number of iterations. Smaller `ϵ`
usually produces better results, but default value should be OK for
all cases. `noise` is an optional array of booleans which contains
initial approximation.

# References
1. A. Cherkasov, A. Ananev, Adaptive phase-retrieval stochastic
   reconstruction with correlation functions: Three-dimensional images
   from two-dimensional cuts, Phys. Rev. E, 104, 3, 2021

See also: [`autocorrelation`](@ref).
"""
function phaserec(target :: AbstractArray{<: AbstractFloat}, arraysize;
                  radius   = 0.6,
                  α        = 0.997,
                  maxsteps = 300,
                  noise  :: Maybe{AbstractArray{Bool}} = nothing,
                  bc     :: BoundaryConditions = Periodic())
    p = porosity(target, arraysize, bc)
    @assert 0 < p < 1

    # Make initial guess for the reconstruction
    guess = isnothing(noise) ? initial_guess(arraysize, p) : noise
    recon  = maybe_pad(guess, bc)

    # Make Gaussian low-pass filter
    filter = make_filter(size(recon), radius)
    # Cost function based on S₂ map
    normfn(corr, rec) = norm(corr - autocorrelation(rec))
    initnorm = normfn(target, recon)
    oldn = 1.0

    for steps in 1:maxsteps
        # Restore S₂ function
        gray = replace_abs(recon, target)

        diff = maximum(abs.(gray)) - minimum(abs.(gray))
        gray = gray + 0.002 * diff * rand(ComplexF64, size(gray))

        # Apply low-pass filter
        gray = filter .* gray
        # Convert to two-phase image
        gray_spatial = maybe_cut_padding(irfft(gray, size(recon, 1)), bc)
        candidate = maybe_pad(threshold(gray_spatial, p), bc)

        n = normfn(target, candidate) / initnorm
        if n < oldn
            # Candidate is better than the previous reconstruction, accept
            recon = candidate
        else
            r = rand()
            thr = exp(-500(n - oldn) / α^steps)
            @info "Rand = $(r), Threshold = $(thr)"
            # Candidate is worse than the previous reconstruction,
            # accept by random choice.
            if r < thr
                recon = candidate
            end
        end
        n = normfn(target, recon) / initnorm

        if (rem(steps, 10) == 1)
            @info "Cost = $(n)"
        end

        oldn = n
    end

    return maybe_cut_padding(recon, bc), normfn(target, recon) / initnorm
end

"""
    phaserec(array; radius = 0.6, maxsteps = 300, ϵ = 10^-5[, noise])

Reconstruct `array` which must be an array of booleans. This is
equivalent to running `phaserec(autocorrelation(array), size(array); ...)`.
"""
phaserec(array :: AbstractArray{Bool};
         radius   = 0.6,
         α        = 0.997,
         maxsteps = 300,
         noise    = nothing,
         bc       = Periodic()) =
    phaserec(autocorrelation(maybe_pad(array, bc)), size(array);
             radius, α, maxsteps, noise, bc)
