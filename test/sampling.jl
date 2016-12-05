# sampling.jl
using Base.Test
using Wavelets
# ============= down- and upsampling ================
println("sampling: downsampling, upsampling and filters ...")


# Compare the efficient implemention of filtering with downsampling to a simple
# brute-force computation.
for M = 1:5
    for L = 1:10
        fslen = Sampling.fstate_length_downsampling(L, M)
        @test fslen ==  ceil(Int, (L-1)/M)
        filter = rand(L)
        for Ny = 100:100+lcm(L,M)
            y1 = zeros(Ny)
            y2 = similar(y1)
            Nx = L+(Ny-1)*M
            x = rand(Nx)
            Sampling.filter_and_downsample!(y1, x, filter, M)
            Sampling.filter_and_downsample_bruteforce!(y2, x, filter, M)
            @test norm(y1-y2) < 1e-10
        end
    end
end

# Compare the efficient implemention of filtering with upsampling to a simple
# brute-force computation.
for M = 1:5
    for L = 1:10
        filter = rand(L)
        for Nx = 100:100+lcm(L,M)
            for M0 = 0:M-1
                x = rand(Nx)
                y1 = zeros(M*Nx-M*L)
                y2 = similar(y1)
                fill!(y2, 0)
                Sampling.upsample_and_filter_bruteforce!(y1, x, filter, M, M0)
                Sampling.upsample_and_filter!(y2, x, filter, M, M0)
                @test norm(y1-y2) < 1e-10
            end
        end
    end
end
