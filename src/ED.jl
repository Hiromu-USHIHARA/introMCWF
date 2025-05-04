module sparse_transform_functions
    # 2Dのindex (m, n)に対し, M=N(m-1)+(n-1)+1を返す
    using SparseArrays
    function _2Dind_to_1Dind(m, n, N) # N=2^nsite
        return N*(m-1)+n
    end
    # 1Dのindex Mに対し, m, nを返す
    function _1Dind_to_2Dind(M, N)
        return (M-1)÷N + 1, (M-1)%N + 1
    end
    # %%
    # 2Dのmatrixに対し, 1Dの配列を返す
    function _2Dmat_to_1Darray(O_mat, N)
        result=spzeros(ComplexF64, N^2)
        for II ∈ 1:N^2
            m, n = _1Dind_to_2Dind(II, N)
            result[II] = O_mat[m, n]
        end
        return result
    end
    # 1Dの配列に対し, 2Dのmatrixを返す
    function _1Darray_to_2Dmat(O_array, N)
        result = spzeros(ComplexF64, N, N)
        for II ∈ 1:N^2
            m, n = _1Dind_to_2Dind(II, N)
            result[m, n] = O_array[II]
        end
        return result
    end
    # %%
    function matrix_w_energy_basis(O_mat, N, vs)
        @assert size(O_mat) == (N, N) "size of O_mat ($(size(O_mat))) is inconsistent with N ($N)"
        result = spzeros(typeof(O_mat[begin, begin]), N, N)
        for n ∈ 1:N
            vn = vs[:, n]
            for m ∈ 1:N
                vm = vs[:, m]
                result[m, n] = vm'*O_mat*vn
            end
        end
        return result
    end
end
# %%
module sparse_Liouville_space
    using LinearAlgebra, KrylovKit, SparseArrays
    using Main.sparse_transform_functions
    ⊗(x, y)=kron(x, y)
    # %%
    function super_commutator(X, Y)
        X⊗transpose(Y) - Y⊗transpose(X)
    end
    #
    function super_anti_commutator(X, Y)
        X⊗transpose(Y) + Y⊗transpose(X)
    end
    #
    function dissipator_in_Liouville_space(γs, ℒs)
        @assert length(γs)==length(ℒs) "γs ($γs) & ℒs ($ℒs) must have same length"
        n = length(γs)
        𝕀 = sparse(I, size(ℒs[1]))
        #
        L = ℒs[1]
        result = γs[1]*(L⊗conj(L) - super_anti_commutator(L' * L, 𝕀)/2.)
        for i ∈ 2:n
            L = ℒs[i]
            result += γs[i]*(L⊗conj(L) - super_anti_commutator(L' * L, 𝕀)/2.)
        end
        return result
    end
    #
    function unitary_in_Liouville_space(H)
        𝕀 = sparse(I, size(H))
        return -im*super_commutator(H, 𝕀)
    end
    #
    function GKSL_ED_in_Liouville_space(H, γs, ℒs)
        𝔘 = unitary_in_Liouville_space(H)
        𝔇 = dissipator_in_Liouville_space(γs, ℒs)
        𝔏 = 𝔘 + 𝔇
        # es, vs = eigen(𝔏)
        es, vs, _ = eigsolve(𝔏, 3, :LR)
        @show es
        return es, vs
    end
    #
    function GKSL_NESS_in_Liouville_space(H, γs, ℒs)
        es, vs = GKSL_ED_in_Liouville_space(H, γs, ℒs)
        # @show inds_NESS = findall(x -> abs(x)<=1.0e-10, es)
        # result = transform_functions._1Darray_to_2Dmat(sum(vs[:, i] for i ∈ inds_NESS), size(H)[1]) # ρ_NESS
        result = sparse_transform_functions._1Darray_to_2Dmat(vs[1], size(H, 1)) # ρ_NESS
    end
    # 
    function time_evol_ED(H, γs, Ls, ρ0, ts)
        ℒ = unitary_in_Liouville_space(H) + dissipator_in_Liouville_space(γs, Ls)
        esL, vsL = eigen(Matrix(ℒ))
        #
        Vinv_ρ0 = vsL \ sparse_transform_functions._2Dmat_to_1Darray(ρ0, size(H, 1))
        @show sum(Vinv_ρ0[i]*vsL[:,i] for i ∈ 1:size(vsL, 2)) ≈ sparse_transform_functions._2Dmat_to_1Darray(ρ0, size(H, 1))
        ρs = [sparse_transform_functions._1Darray_to_2Dmat(sum(exp(esL[i]*t)* Vinv_ρ0[i] * vsL[:,i] for i ∈ 1:(size(H,1)^2)), size(H, 1)) for t ∈ ts]
        return ρs
    end
end

