function isPhase(J; N, J1, nsweeps, Ts)
    data = []
    spin_config = fill([1, 0, 0], N, N) .* ((-1) .^ ((1:10)' .+ (1:10)))

    println("Annealing J2/J1 = $(J)\n--------------")
    @time for T in Ts
        println("T = $(T)")
        push!(data, simulate!(spin_config, T, J1, J .* J1, nsweeps, 0.25))
    end
    println("--------------")
    println("Mππ: $(data[end][3])")

    return !isapprox(data[end][3], 0., atol = 1e-4)
end

function binarySearch(low, high; kwargs...)
    mid = 0
    
    for i in 1:10
        mid = (low + high)/2
        
        if isPhase(mid; kwargs...)
            low = mid
        else
            high = mid
        end
    end

    return mid
end

res = binarySearch(0., 0.6; N = 10, J1 = 4, nsweeps = 10000, Ts = 4.:-0.1:0.1)
println("Jc = $(res)")