function read_data(name, dims, side)
    s = Tuple(side for _ in 1:dims)
    array = zeros(UInt8, s)

    open(name) do io
        read!(io, array)
    end

    return BitArray(array)
end

diff(a1, a2) = norm(two_point(a1) - two_point(a2))

function test_rec(array)
    rec1 = phaserec(array; maxsteps = 1)
    rec2 = phaserec(array; maxsteps = 100)

    @test diff(rec1, array) > diff(rec2, array)
end

@testset "Reconstruct ceramic" begin
    test_rec(read_data("ceramic200", 3, 200))
end

@testset "Reconstruct soil" begin
    test_rec(read_data("soil160", 3, 160))
end
