# Exhaustive search for factorization of elements of SL4(GF(2)) with upper right block in SL2 using block unitriangular matrices 

using SimpleGF2, LinearAlgebra, StaticArrays

function SL4GF2()
    T = GF2 # type
    SM = SMatrix{4, 4, T, 16}

    M = MMatrix{4,4,T}(undef)

    SL = Vector{SM}() # SL4(GF(2)) with upper right block in SL2(GF(2))
    D5 = Vector{SM}() # five-fold products of block unitriangular matrices

    M2 = [ @SMatrix[a1 a2; a3 a4] for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1 ] # M2(GF(2))

    BL4 = SM[ SMatrix{4,4}(T.([1 0 0 0; 0 1 0 0; a1 a2 1 0; a3 a4 0 1])) for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1 ] # block lower unitriangular
    BU4 = SM[ SMatrix{4,4}(T.([1 0 a1 a2; 0 1 a3 a4; 0 0 1 0; 0 0 0 1])) for a1 in 0:1, a2 in 0:1, a3 in 0:1, a4 in 0:1 ] # block upper unitriangular

    for a ∈ 1:16
        M[1:2,1:2] = M2[a]
        for b ∈ 1:16
            AB = BL4[a] * BU4[b] # two-fold products BL4 x BU4
            M[1:2,3:4] = M2[b]
            for c ∈ 1:16
                ABC = AB * BL4[c] # three-fold products BL4 x BU4 x BL4
                M[3:4,1:2] = M2[c]
                for d ∈ 1:16
                    ABCD = ABC * BU4[d] # four-fold products BL4 x BU4 x BL4 x BU4
                    M[3:4,3:4] = M2[d]
                    let M = SMatrix(M) # done with modifying it
                        if one(T) == det(M) == det(SMatrix{2,2}(view(M, 1:2,3:4)))
                            push!(SL, M) # if matrix in SL4(GF(2)) with upper right block in SL2(GF(2)), store in SL
                        end
                    end
                    for e ∈ 1:16
                        ABCDE = ABCD * BL4[e] # five-fold products BL4 x BU4 x BL4 x BU4 x BL4
                        push!(D5, ABCDE)  # store five-fold products in D5
                    end
                end
            end
        end
    end

    # D5 = unique(D5) # remove copies from D5
    setdiff(SL,D5) # output all matrices in SL4(GF(2)) with upper right block in SL2(GF(2)) that are not a five-fold product BL4 x BU4 x BL4 x BU4 x BL4
    # Output is []. QED
end
