import JLD: save

N = 10
J1 = 4

# Js = [0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]
Js = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05]

nsweeps = 5000
Ts = 4.:-0.1:0.

@time for J in Js
    data = []
    spin_config = fill([1, 0, 0], N, N) .* ((-1) .^ ((1:10)' .+ (1:10)))

    println("Annealing J2/J1 = $(J)\n--------------")
    for T in Ts
        println("T = $(T)")
        push!(data, simulate!(spin_config, T, J1, J .* J1, nsweeps, 0.25))
    end
    println("--------------")
    save("./data/ising_J2J1_$(J)_N_$(N)_sweeps_$(nsweeps)_AFM_init.jld", "Ts", Ts, "J", J, "N", N, "data", data) 
end

