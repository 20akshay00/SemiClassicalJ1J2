import Jackknife: estimate
import LinearAlgebra: norm
import StatsBase: mean 

random_spin() = [2 * Int(rand() < 0.5) - 1, 0, 0]

function is_spin(spin_config, site, N)
    s = spin_config[mod1(site[1], N), mod1(site[2], N)]
    return length(s) == 3
end

energy_singlet(J, beta) = ((3*J)/4) * ((exp(-beta * J) - 1)/(3 * exp(-beta * J) + 1))
dot(v1, v2) = (length(v1) == 3 && length(v2) == 3) ? v1' * v2 : 0

function dist(site1, site2, N)
    rx, ry = abs.(site1 .- site2)
    dist_x = (rx >= N/2) ? (N - rx) : rx
    dist_y = (ry >= N/2) ? (N - ry) : ry

    if iszero(dist_x) || iszero(dist_y)
        dist_tot = 1
    elseif isone(dist_x) && isone(dist_y)
        dist_tot = 2
    else
        print("Oh no!")
    end

    return dist_tot
end

function energy_spin(spin_config, J1, J2)
    N = size(spin_config, 1)
    energy = 0 
    for i in 1:N
        for j in 1:N
            energy += (J1/4) * (dot(spin_config[i, j], spin_config[mod1(i + 1, N), j]) + dot(spin_config[i, j], spin_config[i, mod1(j + 1, N)]))
            energy += (J2/4) * (dot(spin_config[i, j], spin_config[mod1(i + 1, N), mod1(j + 1, N)]) + dot(spin_config[i, j], spin_config[mod1(i + 1, N), mod1(j - 1, N)]))
        end
    end

    return energy
end

function get_energy(spin_config, J1, J2, beta)
    N = size(spin_config, 1)
    nJ1, nJ2 = 0, 0

    for i in 1:N
        for j in 1:N
            if(length(spin_config[i, j]) == 2)
                d = dist([i, j], spin_config[i, j], N)

                nJ1 += (d == 1)
                nJ2 += (d == 2)
            end
        end
    end

    return energy_spin(spin_config, J1, J2) + nJ1/2 * energy_singlet(J1, beta) + nJ2/2 * energy_singlet(J2, beta)
end

function make_dimer!(site, spin_config, J1, J2, beta)
    N = size(spin_config, 1)
    unpaired_list = []
    ix, iy = site 

    for nbr in [[ix+1, iy], [ix-1, iy], [ix, iy+1], [ix, iy-1], [ix-1, iy-1], [ix-1, iy+1], [ix+1, iy-1], [ix+1, iy+1]]
        if is_spin(spin_config, site, N)
            push!(unpaired_list, nbr)
        end
    end
        
    if !iszero(length(unpaired_list)) 
        unpaired_site = mod1.(rand(unpaired_list), N)
        jx, jy = unpaired_site 

        d = dist(site, unpaired_site, N)
        J = (d == 1) ? J1 : J2

        E1 = energy_singlet(J, beta)
        E2 = - (J1/4) * (dot(spin_config[ix, iy], spin_config[mod1(ix - 1, N), iy]) + dot(spin_config[ix, iy], spin_config[mod1(ix + 1, N), iy]) + dot(spin_config[ix, iy], spin_config[ix, mod1(iy - 1, N)]) + dot(spin_config[ix, iy], spin_config[ix, mod1(iy + 1, N)]) + dot(spin_config[jx, jy], spin_config[mod1(jx - 1, N), jy]) + dot(spin_config[jx, jy], spin_config[mod1(jx + 1, N), jy]) + dot(spin_config[jx, jy], spin_config[jx, mod1(jy + 1, N)]) + dot(spin_config[jx, jy], spin_config[jx, mod1(jy - 1, N)]))
        E3 = - (J2/4) * (dot(spin_config[ix, iy], spin_config[mod1(ix - 1, N), mod1(iy - 1, N)]) + dot(spin_config[ix, iy], spin_config[mod1(ix + 1, N), mod1(iy - 1, N)]) + dot(spin_config[ix, iy], spin_config[mod1(ix - 1, N), mod1(iy + 1, N)]) + dot(spin_config[ix, iy], spin_config[mod1(ix + 1, N), mod1(iy + 1, N)]) + dot(spin_config[jx, jy], spin_config[mod1(jx - 1, N), mod1(jy - 1, N)]) + dot(spin_config[jx, jy], spin_config[mod1(jx + 1, N), mod1(jy - 1, N)]) + dot(spin_config[jx, jy], spin_config[mod1(jx - 1, N), mod1(jy + 1, N)]) + dot(spin_config[jx, jy], spin_config[mod1(jx + 1, N), mod1(jy + 1, N)]))
        E4 = (J/4) * (dot(spin_config[ix, iy], spin_config[jx, jy]))
        dE = E1 + E2 + E3 + E4
        
        if rand() <= exp(-beta * dE)
            spin_config[ix, iy] = [jx, jy]
            spin_config[jx, jy] = [ix, iy]
        else 
            dE = 0
        end
    else 
        dE = 0 
    end

    return dE    
end

function break_dimer!(site, spin_config, J1, J2, beta)
    N = size(spin_config, 1)
    d = dist(site, spin_config[site...], N)
    J = (d == 1) ? J1 : J2

    spin_i, spin_j = random_spin(), random_spin()
    ix, iy = site
    jx, jy = spin_config[site...]

    E1 = - energy_singlet(J, beta)
    E2 = (J1/4) * (dot(spin_i, spin_config[mod1(ix - 1, N), iy]) + dot(spin_i, spin_config[mod1(ix + 1, N), iy]) + dot(spin_i, spin_config[ix, mod1(iy - 1, N)]) + dot(spin_i, spin_config[ix, mod1(iy + 1, N)]) + dot(spin_j, spin_config[mod1(jx - 1, N), jy]) + dot(spin_j, spin_config[mod1(jx + 1, N), jy]) + dot(spin_j, spin_config[jx, mod1(jy + 1, N)]) + dot(spin_j, spin_config[jx, mod1(jy - 1, N)]))
    E3 = (J2/4) * (dot(spin_i, spin_config[mod1(ix - 1, N), mod1(iy - 1, N)]) + dot(spin_i, spin_config[mod1(ix + 1, N), mod1(iy - 1, N)]) + dot(spin_i, spin_config[mod1(ix - 1, N), mod1(iy + 1, N)]) + dot(spin_i, spin_config[mod1(ix + 1, N), mod1(iy + 1, N)]) + dot(spin_j, spin_config[mod1(jx - 1, N), mod1(jy - 1, N)]) + dot(spin_j, spin_config[mod1(jx + 1, N), mod1(jy - 1, N)]) + dot(spin_j, spin_config[mod1(jx - 1, N), mod1(jy + 1, N)]) + dot(spin_j, spin_config[mod1(jx + 1, N), mod1(jy + 1, N)]))
    E4 = - (J/4) * (dot(spin_i, spin_j))
    dE = E1 + E2 + E3 + E4
    
    if (rand() <= exp(-beta * dE))
        spin_config[ix, iy] = spin_i
        spin_config[jx, jy] = spin_j
    else
        dE = 0
    end

    return dE
end

function rotate_spin!(site, spin_config, J1, J2, beta)
    N = size(spin_config, 1)
    spin_new = random_spin()
    ix, iy = site 

    dS = spin_new .- spin_config[site...]

    dE1= (J1/4) * (dot(dS, spin_config[mod1(ix - 1, N), iy]) + dot(dS, spin_config[mod1(ix + 1, N), iy]) + dot(dS, spin_config[ix, mod1(iy - 1, N)]) + dot(dS, spin_config[ix, mod1(iy + 1, N)])) 
    dE2= (J2/4) * (dot(dS, spin_config[mod1(ix - 1, N), mod1(iy - 1, N)]) + dot(dS, spin_config[mod1(ix + 1, N), mod1(iy - 1, N)]) + dot(dS, spin_config[mod1(ix - 1, N), mod1(iy + 1, N)]) + dot(dS, spin_config[mod1(ix + 1, N), mod1(iy + 1, N)]))
    
    dE = dE1 + dE2

    if rand() < exp(-beta * dE)
        spin_config[site...] = spin_new
    else
        dE = 0
    end

    return dE
end

function sweep!(spin_config, J1, J2, beta, f, E)
    N = size(spin_config, 1)

    for i in 1:N^2
        site = rand(1:N, 2)

        if is_spin(spin_config, site, N)
            if rand() < f
                E += make_dimer!(site, spin_config, J1, J2, beta)
            else
                E += rotate_spin!(site, spin_config, J1, J2, beta)
            end
        else
            E += break_dimer!(site, spin_config, J1, J2, beta)
        end
    end

    free_spins = []
    n_dimer = 0 

    for i in 1:N
        for j in 1:N
            if is_spin(spin_config, [i, j], N)
                push!(free_spins, [i, j])
            else
                n_dimer += 1
            end
        end
    end

    S, magππ, mag0π, magπ0 = zeros(3), zeros(3), zeros(3), zeros(3)

    if !iszero(length(free_spins))
        for site in free_spins
            spin = spin_config[site...]

            S += spin
            magππ += spin * (-1)^(sum(site))
            magπ0 += spin * (-1)^(site[1])
            mag0π += spin * (-1)^(site[2])
        end
    end

    return E, S, magππ, mag0π, magπ0, n_dimer
end

function simulate!(spin_config, T, J1, J2, nsweeps, f)
    N = size(spin_config, 1)
    beta = 1/T
    E0 = get_energy(spin_config, J1, J2, beta)

    # observables
    neq_sweeps = nsweeps ÷ 10
    nsweeps = nsweeps - neq_sweeps

    E = Vector{Float64}(undef, nsweeps)
    S = Vector{Vector{Int64}}(undef, nsweeps)
    Mππ = Vector{Vector{Int64}}(undef, nsweeps)
    M0π = Vector{Vector{Int64}}(undef, nsweeps)
    Mπ0 = Vector{Vector{Int64}}(undef, nsweeps)
    nD = Vector{Int64}(undef, nsweeps)

    # equilibriation
    for i in 1:neq_sweeps
        E0, _, _, _, _, _ = sweep!(spin_config, J1, J2, beta, f, E0)
    end

    # first sweep to pass the energy
    E[1], S[1], Mππ[1], M0π[1], Mπ0[1], nD[1] = sweep!(spin_config, J1, J2, beta, f, E0)

    # gather data 
    for i in 2:nsweeps
        E[i], S[i], Mππ[i], M0π[i], Mπ0[i], nD[i] = sweep!(spin_config, J1, J2, beta, f, E[end])
    end
        
    # store data
    return estimate.([mean], [E, norm.(S), norm.(Mππ), norm.(M0π), norm.(Mπ0), nD])
end