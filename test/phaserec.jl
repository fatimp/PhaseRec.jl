function read_data(name, dims, side)
    s = Tuple(side for _ in 1:dims)
    array = zeros(UInt8, s)

    open(name) do io
        read!(io, array)
    end

    return BitArray(array)
end

function test_rec(array)
    _, imp1 = phaserec(array; maxsteps = 1)
    _, imp2 = phaserec(array; maxsteps = 100)

    @test imp2 < imp1
end

@testset "Reconstruct ceramic" begin
    test_rec(read_data("ceramic200", 3, 200))
end

@testset "Reconstruct soil" begin
    test_rec(read_data("soil160", 3, 160))
end
