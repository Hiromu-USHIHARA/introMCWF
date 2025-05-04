module trajectory_approach
    using LinearAlgebra, SparseArrays, Random
    # %%
    function time_evol_trajectory(ψ0, dt, tmax, H, Ls)
        Heff = H - 0.5*im*sum(L'*L for L in Ls)
        # es, vs = eigen(Matrix(Heff))
        ts = 0:dt:tmax
        ψs = [Vector{ComplexF64}(ψ0)]
        for i ∈ 1:(length(ts)-1)
            ψ = ψs[i]
            # Vinvψ =  vs\ψ
            # ψ = sum(exp(-im*dt*es[i])*Vinvψ[i]*vs[:,i] for i ∈ 1:size(vs,2))
            ψ = (I(size(H,1)) - im*dt*Heff)*ψ
            δps = [dt*real(ψs[i]'*Ls[j]'*Ls[j]*ψs[i]) for j ∈ 1:length(Ls)]
            δp = sum(δps)
            r1 = rand()
            if r1 > δp
                # ψ = ψ/√(1-δp)
                ψ = ψ/√(ψ'*ψ)
            else
                r2 = rand()
                j = findfirst(cumsum(δps/δp) .> r2)
                ψ = Ls[j]*ψs[i]
                ψ = ψ/√(ψ'*ψ)
            end
            append!(ψs, [ψ])
        end
        return ts, ψs
    end  
end
