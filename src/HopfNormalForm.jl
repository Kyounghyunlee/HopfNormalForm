"""
 Centermanifold computation and simplest normal form of Hopf bifurcation.
 This code is available to compute centermanifold and normal form of only semisimple Jacobian case (i.e. diagonalizable Jacobian).
 For nonsemisimple case, algorithm should be modified.

 SymEngine package is used.

 For theoretical background of Centermanifold computation see:
    http://www.orcca.on.ca/TechReports/TechReports/2000/TR-00-15.pdf
 for simplest normal form of Hopf bifurcation see:
    https://iopscience.iop.org/article/10.1088/0951-7715/16/1/317
"""
module HopfNormalForm

using SymEngine, OffsetArrays

export calculate_normal_form, verify

################### Basic functions for symbolic computation in Julia ###################
function coeff2(f::Basic, e::Array{Basic,1}, N)::Basic #coefficient of multivariate Polynomial for univariate case just use diff
    # f is function of plinomials of variable e and N is index vector of power series. for example,
    # if f=e(1)^j*e(2)^k -> I=[j,k]
    v = f
    for i in eachindex(N)
        if N[i] == 0
            v *= 1
        else
            v = diff(v, e[i], N[i]) / factorial(N[i])
        end
    end
    for j in eachindex(N)
        v = subs(v, e[j] => 0)
    end
    return v
end

function recursive_apply(f, x::Basic)::Basic
    cls = SymEngine.get_symengine_class(x)
    args = [recursive_apply(f, arg) for arg in get_args(x)]
    if cls == :Add
        return sum(args)
    elseif cls == :Mul
        return prod(args)
    elseif cls == :Pow
        return args[1]^args[2]
    elseif cls in [:RealDouble, :RealMPFR, :ComplexDouble, :ComplexMPFR]
        return f(x)
    else
        return x
    end
end

function chop_expr(x::Basic, zero_tol::Real)::Basic #chopping very small coefficients to simplify the equaion expression
    return recursive_apply(x -> (abs(N(x)) <= zero_tol) ? zero(typeof(x)) : x, x)
end

function recursive_real(x::Basic)::Basic # Real part of complex symbols. All symbol variables are assumed to be real.
    # To define a complex symbol variable c, express it as c=a+b*I.
    return recursive_apply(real, x)
end

function recursive_imag(x::Basic)::Basic # Imaginary part of complex symbols. All symbol variables are assumed to be real.
    # To define a complex symbol variable c, express it as c=a+b*I.
    return recursive_apply(imag, x)
end

function verify(
    f::Array{Basic},
    T::Array{Basic},
    v::Array{Basic},
    norder::Integer,
    cm::Integer,
    m3::Integer,
    n::Integer,
)::Array{Basic}
    # verify the result of Centermanifold computation upto norder
    # f: function input, T: computed centermanifold, v: state-space variables, norder: oreder of computation,
    # cm: dimension of centermanifold, m3: dimension of stable-unstable manifold
    G2 = Array{Basic}(undef, m3)
    dhf = Array{Basic}(undef, m3)
    dh = Array{Basic}(undef, m3)
    F = Array{Basic}(undef, cm)
    veri = Array{Basic}(undef, m3)
    veri_o = Array{Basic}(undef, m3, norder - 1)
    @vars ϵ
    # Compute F
    for i = 1:cm
        F[i] = f[i] - Eign[i] * v[i]
        for m = 1:n
            if m > cm
                F[i] = subs(F[i], v[m] => T[m-cm])
                F[i] = expand(F[i])
                F[i] = chop_expr(F[i], 1e-10)
            end
        end
    end
    # Compute dh*f
    for i = 1:m3
        dhf[i] = 0
        for l = 1:cm
            dhf[i] += diff(T[i], v[l]) * F[l]
        end
    end
    # Compute g-dh*f
    for k = 1:m3
        G[k] = f[cm+k] - Eign[cm+k] * v[cm+k]
        for m = 1:n
            if m > cm
                G[k] = subs(G[k], v[m] => T[m-cm])
            end
        end
        G2[k] = G[k] - dhf[k]
    end
    # compute dh(x)J_cx-J_s*h(x)
    for i = 1:m3
        dh[i] = 0
        for l = 1:cm
            dh[i] += diff(T[i], v[l]) * Eign[l] * v[l]
        end
        dh[i] = dh[i] - Eign[i+cm] * T[i]
    end
    # veri[i]= g-dh*f - {dh(x)J_cx-J_s*h(x)}
    # We have to check order j=2,3,…,norder of veri[i] is zero
    for i = 1:m3
        veri[i] = dh[i] - G2[i]
        veri[i] = expand(veri[i])
        veri[i] = chop_expr(veri[i], 1e-10)
    end
    # returning value veri_o[i,j] is j-th order of veri[i]
    for j = 2:norder
        for i = 1:m3
            for m = 1:cm
                veri[i] = subs(veri[i], v[m] => ϵ * v[m])
            end
            veri_o[i, j-1] = subs(diff(veri[i], ϵ, j), ϵ => 0) / factorial(j)
        end
    end
    return veri_o
end

#==
Note that complex conjugate pair of eigenvalue corresponding to center subspace is transfromed to +-i.
This is equivalent to rescailing time.

################### Define the index for recursive computation ###################
To compute centermanifold of higher dimension n>3, you should modify this part.
index of coefficient arrays are [k,l,m,n] for centermanifold of dimension 3
k is for k-th part of function g=[g_1,g_2,…,g_k,…,g_m3], l,m,n is for coefficient of polynomial of v_1^l* v_2^m* v_3^n
of power series of g_k
==#

function calculate_normal_form(
    f::Vector{Basic},
    v::Vector{Basic},
    norder::Integer,
    cm::Integer,
    m1::Integer,
    m2::Integer,
    m3::Integer,
    n::Integer,
    ω_h::AbstractFloat,
    zero_tol::AbstractFloat
)
    G = Array{Basic}(undef, m3);
    T = Array{Basic}(undef, m3);
    d = Array{Basic}(undef, m3);
    F = Array{Basic}(undef, cm);

    a = Array{Basic}(undef, cm)
    RD = Array{Basic}(undef, cm, norder); # Reduced dynamics on the center manifold by order
    RD_sum = Array{Basic}(undef, cm); # Reduced dynamics on the center manifold summation of all orders

    @vars ϵ
    # ϵ is symbolic variable which is used in the algorithm which does not have any physical meaning

    if cm == 2 # set multidimensional array for saving coefficints of the polynomial
    # for cm>=3 you should modify this part
        B = OffsetArray{Basic}(undef, 1:m3, 0:norder, 0:norder)
        H = OffsetArray{Basic}(undef, 1:m3, 0:norder, 0:norder)
    elseif cm == 3
        B = OffsetArray{Basic}(undef, 1:m3, 0:norder, 0:norder, 0:norder)
        H = OffsetArray{Basic}(undef, 1:m3, 0:norder, 0:norder, 0:norder)
        A = OffsetArray{Basic}(undef, 1:2, 0:norder, 0:norder, 0:norder)
    end

    Eign = Array{Basic}(undef, n); # Array to save eigenvalues

    ########## Collect Eigenvalues ############
    # Eigenvalues are collected from the input function
    for i = 1:n
        Eign[i] = diff(f[i], v[i])
        for k = 1:n
            Eign[i] = subs(Eign[i], v[k] => 0)
        end
    end

    ########## Initializing the centermanifold and function d (i.e. Dh(x)f(x,h(x))) ######
    for k = 1:m3
        T[k] = 0
        d[k] = 0
    end

    ############ Compute g(x,h(x)) ##############
    for j = 2:norder
        for k = 1:m3
            G[k] = f[cm+k] - Eign[cm+k] * v[cm+k]
            for m = 1:n
                if m > cm
                    G[k] = subs(G[k], v[m] => T[m-cm])
                end
            end
            G[k] = G[k] - d[k]
            for m = 1:cm
                G[k] = subs(G[k], v[m] => ϵ * v[m])
            end
            G[k] = subs(diff(G[k], ϵ, j), ϵ => 0) / factorial(j) # Collect the j-th order of the G[k]
            G[k] = expand(G[k])
            G[k] = chop_expr(G[k], zero_tol)

    ######## Define the index of polynomial ########
            # ind_list = Vector{Vector{Int}}()
            # if cm == 2
            #     for l = 0:j
            #         m = j - l
            #         push!(ind_list, [l, m])
            #     end
            # elseif cm == 3
            #     for l = 0:j, m = 0:j-l
            #         nn = j - l - m
            #         push!(ind_list, [l, m, nn])
            #     end
            # end
            # ind_l = Vector{Vector{Int}}()
            # for p = 1:length(ind_list)
            #     push!(ind_l, [k; ind_list[p]])
            # end
            if cm == 2
                ind_l = ((k, l, j - l) for l = 0:j)
            elseif cm == 3
                ind_l = ((k, l, m, j - l - m) for l = 0:j for m = 0:j-l)
            end
            for arg in ind_l
                xs = one(Basic)
                δ = zero(Basic)
                B[CartesianIndex(arg)] = coeff2(G[arg[1]], v, arg[2:end])
                for i = 1:cm
                    B[CartesianIndex(arg)] =
                        subs(B[CartesianIndex(arg)], v[i] => 0)
                end
                for i = 1:cm
                    δ += Eign[i] * arg[1+i]
                    xs *= v[i]^arg[1+i]
                end
                δ = δ - Eign[cm+arg[1]]
                H[CartesianIndex(arg)] = B[CartesianIndex(arg)] / δ
                H[CartesianIndex(arg)] = expand(H[CartesianIndex(arg)])
                H[CartesianIndex(arg)] =
                    chop_expr(H[CartesianIndex(arg)], zero_tol)
                T[arg[1]] += H[CartesianIndex(arg)] * xs
                T[arg[1]] = expand(T[arg[1]])
            end
        end
    ############ Compute Dh(x)f(x,h(x)) ##############
        for i = 1:cm
            F[i] = f[i] - Eign[i] * v[i]
            for m = 1:n
                if m > cm
                    F[i] = subs(F[i], v[m] => T[m-cm])
                    F[i] = expand(F[i])
                    F[i] = chop_expr(F[i], zero_tol)
                end
            end
            d[i] = 0
        end
        for i = 1:m3, l = 1:cm
            d[i] += diff(T[i], v[l]) * F[l]
            d[i] = expand(d[i])
            d[i] = chop_expr(d[i], zero_tol)
        end
    end
    ########## Computing reduced dynamics on the centermanifold ##########
    for k = 1:cm, j = 1:norder
        if j == 1
            RD[k, j] = Eign[k] * v[k]
        else
            RD[k, j] = f[k] - Eign[k] * v[k]
            for m = 1:n
                if m > cm
                    RD[k, j] = subs(RD[k, j], v[m] => T[m-cm])
                end
            end
            for m = 1:cm
                RD[k, j] = subs(RD[k, j], v[m] => ϵ * v[m])
            end
            RD[k, j] = subs(diff(RD[k, j], ϵ, j), ϵ => 0) / factorial(j) # Collect the j-th order of Reduced dynmics
            RD[k, j] = expand(RD[k, j])
            RD[k, j] = chop_expr(RD[k, j], zero_tol)
        end
    end
    for n = 1:cm
        RD_sum[n] = sum(RD[n, :])
    end
    ########## Simplest Normal form computation of Hopf bifurcation ##########

    A[1, 1, 0, 1] = recursive_real(coeff2(RD_sum[2], v, [1, 1, 0]));
    A[1, 2, 1, 0] = recursive_real(coeff2(RD_sum[2], v, [0, 2, 1]));
    A[1, 2, 0, 0] = recursive_real(coeff2(RD_sum[2], v, [0, 2, 0]));
    A[2, 1, 1, 0] = recursive_imag(coeff2(RD_sum[2], v, [0, 1, 1]));
    A[2, 2, 0, 0] = recursive_imag(coeff2(RD_sum[2], v, [0, 2, 0]));
    A[1, 1, 1, 0] = recursive_real(coeff2(RD_sum[2], v, [0, 1, 1]));

    S = A[1, 2, 1, 0] - A[1, 2, 0, 0] * A[2, 1, 1, 0] - A[2, 2, 0, 0] * A[1, 1, 1, 0];
    S = expand(S);

    return -A[1, 1, 0, 1] / S
end

# Include example systems
include("Systems.jl")

end # module
