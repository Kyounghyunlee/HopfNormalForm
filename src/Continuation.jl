module Continuation

using DifferentialEquations
using LinearAlgebra
using NLsolve
using ForwardDiff

BLAS.set_num_threads(1) # This prevents stack overfolowerror

function fourier_diff(N::Integer)
    # For a MATLAB equivalent see http://appliedmaths.sun.ac.za/~weideman/research/differ.html
    h = 2π / N
    D = zeros(N, N)
    # First column
    n = ceil(Int, (N - 1) / 2)
    if (N % 2) == 0
        for i in 1:n
            D[i+1, 1] = -((i % 2) - 0.5) * cot(i * h / 2)
        end
    else
        for i in 1:n
            D[i+1, 1] = -((i % 2) - 0.5) * csc(i * h / 2)
        end
    end
    for i in n+2:N
        D[i, 1] = -D[N-i+2, 1]
    end
    # Other columns (circulant matrix)
    for j in 2:N
        D[1, j] = D[N, j-1]
        for i in 2:N
            D[i, j] = D[i-1, j-1]
        end
    end
    return D
end

function periodic_zero_problem(rhs, u, p, ee)
    # Assume that u has dimensions nN + 1 where n is the dimension of the ODE
    n = dimension(rhs)
    N = (length(u) - 1) ÷ n
    D = fourier_diff(N)
    T = u[end] / (2π)  # the differentiation matrix is defined on [0, 2π]
    # Evaluate the right-hand side of the ODE
    res = Vector{Float64}(undef, n * N + 1)
    for i in 1:n:n*N
        res[i:i+n-1] = T * rhs(@view(u[i:i+n-1]), p, 0)  # use a view for speed
    end
    # Evaluate the derivative; for-loop for speed equivalent to u*Dᵀ
    for (i, ii) in pairs(1:n:n*N)
        # Evaluate the derivative at the i-th grid point
        for (j, jj) in pairs(1:n:n*N)
            for k = 1:n
                res[ii+k-1] -= D[i, j] * u[jj+k-1]
            end
        end
    end
    res[end] = (u[1] - ee)  # phase condition - assumes that the limit cycle passes through the Poincare section u₁=1e-4; using a non-zero value prevents convergence to the trivial solution
    return res
end

function periodic_zero_problem2(rhs, u, p, ee, dv, pu, ds) # Periodic zero problem including Pseudo arc-length equation
    # Assume that u has dimensions nN + 1 where n is the dimension of the ODE
    n = dimension(rhs)
    N = (length(u) - 1) ÷ n
    res = Vector{Float64}(undef, n * N + 2)
    res1 = periodic_zero_problem(rhs, u, p, ee)
    p1 = p[1]
    uu = [u; p1]
    du = uu - pu

    res[1:end-1] = res1
    arclength = norm(transpose(dv) * du)
    res[end] = arclength - ds  #Pseudo-arclength equation
    return res
end

function periodic_zero_J(jeq, em, u, p, ee) # Jacobian of periodic zero problem
    n = dimension(em)
    N = (length(u) - 1) ÷ n
    T = u[end] / (2π)
    J = zeros(n * N + 1, n * N + 1)
    D = fourier_diff(N)
    for i in 1:n:n*N
        J[i:i+n-1, i:i+n-1] = T * jeq(@view(u[i:i+n-1]), p, 0)
        J[i:i+n-1, n*N+1] = em(@view(u[i:i+n-1]), p, 0)
    end
    for (i, ii) in pairs(1:n:n*N)
        for (j, jj) in pairs(1:n:n*N)
            for k = 1:n
                J[ii+k-1, jj+k-1] -= D[i, j]
            end
        end
    end
    J[n*N+1, 1:n*N+1] = zeros(1, n * N + 1)
    J[n*N+1, 1] = 1
    return J
end

function periodic_zero_J2(jeq, jeq2, em, u, p, ee, dv)  # Jacobian of periodic zero problem 2
    n = dimension(em)
    N = (length(u) - 1) ÷ n
    J = zeros(n * N + 2, n * N + 2)
    T = u[end] / (2π)
    J[1:n*N+1, 1:n*N+1] = periodic_zero_J(jeq, em, u, p, ee)
    for i in 1:n:n*N
        J[i:i+n-1, n*N+2] = T * jeq2(@view(u[i:i+n-1]), p, 0)
    end
    J[n*N+2, :] = dv
    return J
end

function LCO_solve(u1, ini_T, eq, jeq, p, ee, z_tol) # Solve the periodic zero problem with Newton methods
    u0 = vcat(u1, ini_T)
    err = []
    for i in 1:40
        res = periodic_zero_problem(eq, u0, p, ee)
        u1 = u0 - periodic_zero_J(jeq, eq, u0, p, ee) \ res
        u0 = u1
        er = transpose(res) * res
        err = vcat(err, er)
        if er < z_tol
            break
        end
    end
    return (u=u0, err=err)
end

function LCO_solve2(u1, ini_T, eq, jeq, jeq2, p, ee, z_tol, dv, pu, ds) # Solve the periodic zero problem 2 with Newton methods
    p0 = p[1]
    u0 = [u1; ini_T]
    uu0 = [u0; p0]
    err = []
    np = p
    for i in 1:40
        res = periodic_zero_problem2(eq, u0, np, ee, dv, pu, ds)
        uu1 = uu0 - periodic_zero_J2(jeq, jeq2, eq, u0, np, ee, dv) \ res
        uu0 = uu1
        u0 = uu0[1:end-1]
        p0 = uu0[end]
        np = vcat(p0, p[2:end])
        er = norm(transpose(res) * res)
        err = vcat(err, er)
        if er < z_tol
            break
        end
    end
    return (u=uu0, err=err)
end

function get_sol(u0, N, ind1, ind2) # convert solution of collocation -> Array form
    # returns phase angle, amplitude
    dim = Int((length(u0) - 1) / N)
    u = u0[1:end-1]
    T = u0[end]
    uu = Array{Float64}(undef, dim, N)
    theta = Array{Float64}(undef, N)
    r = Array{Float64}(undef, N)
    for i in 1:dim
        uu[i, :] = u[i:dim:end]
    end
    for i in 1:N
        theta[i] = atan(uu[ind2, i], uu[ind1, i])
        r[i] = sqrt(uu[ind1, i]^2 + uu[ind2, i]^2)
    end
    return (u=uu, T=T, t=theta, r=r)
end

function get_sol(u0, dim, N, ind1, ind2)  # Convert solution of periodic zero-problem -> array solution and amplitude, phase angle of measured coordinate (z_1,z_2)
    u = u0[1:end-1]
    T = u0[end]
    uu = Array{Float64}(undef, N, dim)
    theta = Array{Float64}(undef, N)
    r = Array{Float64}(undef, N)
    for i in 1:dim
        uu[:, i] = u[i:dim:end]
    end
    for i in 1:N
        theta[i] = atan(uu[i, ind2], uu[i, ind1])
        r[i] = sqrt(uu[i, ind1]^2 + uu[i, ind2]^2)
    end
    return (u=uu, T=T, t=theta, r=r)
end

function get_sol2(sol, ind) # Convert solution of ODE solver to array solution
    lu = length(sol.u)
    u = Vector{Float64}(undef, lu)
    for i in 1:lu
        uv = sol.u[i]
        u[i] = uv[ind]
    end
    return u
end

function MonoD(eq, dim, u, p, t) # Variational equation for monodromy matrix computation (ForwardDiff is used for jacobian computation)
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    n = dim
    v = u[1:Int(n * n)]
    uu = u[Int(n * n + 1):end]
    M = reshape(v, n, n)
    J = jeq(uu, p, t)
    dM = J * M
    dv = vec(dM)
    duu = eq(uu, p, t)
    duu = vec(duu)
    du = [dv; duu]
    return du
end

function Monodromy_compute(eq, u, p, N, ind) # Monodromy matrix computation using numerical integration
    n = dimension(eq)
    eye = 1.0 * Matrix(I, n, n)
    v0 = vec(eye)
    g = get_sol(u, n, N, ind[1], ind[2])
    T = g.T
    uu = g.u
    M = eye
    tl2 = T / N
    for i in 1:N
        uu0 = uu[i, :]
        u0 = [v0; uu0]
        prob = ODEProblem((u, p, t) -> MonoD(eq, n, u, p, t), u0, (0, tl2), p)
        sol = solve(prob, Tsit5(), reltol=1e-14, abstol=1e-14)
        w = sol.u[end]
        m1 = w[1:Int(n * n)]
        M1 = reshape(m1, n, n)
        M = M1 * M
    end
    Eig1 = eigen(M)
    μ = Eig1.values
    return μ
end

function periodic_zero_J_M(jeq, em, u, p) # Jacobian of periodic zero problem
    n = dimension(em)
    N = (length(u) - 1) ÷ n
    T = u[end] / (2π)
    J = zeros(n * N, n * N )
    D = fourier_diff(N)
    
    for i in 1:n:n*N
        J[i:i+n-1, i:i+n-1] = T*jeq(@view(u[i:i+n-1]), p, 0)
    end
    for (i, ii) in pairs(1:n:n*N)
        for (j, jj) in pairs(1:n:n*N)
            for k = 1:n
                J[ii+k-1, jj+k-1] -= D[i, j]
            end
        end
    end
    return J
end
#=
function Monodromy_compute(eq, u, p, N, ind) # Monodromy matrix computation using numerical integration
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    M=periodic_zero_J_M(jeq, eq, u, p)
    Eig1 = eigen(M)
    σ = Eig1.values # Floquet exponent
    μ = exp.(σ* 2π)
    return μ
end
=#
function amp_LCO(U, ind, N) #Extract amplitude from collocation points
    l1 = length(U)
    amp = Vector{Float64}(undef, l1)
    for i in 1:l1
        g = get_sol(U[i], 6, N, 1, 3)
        uu = g.u[:, ind]
        amp[i] = maximum(uu) - minimum(uu)
    end
    return amp
end

function Hopf_point(eq, p) # Compute the Hopf point
    t = 0.0
    u = zeros(dimension(eq))
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    s = nlsolve(p0 -> maximum(real(eigen(jeq(u, p0, t)).values)) - 1e-11, p)
    return s.zero[1]
end

function Hopf_point2(eq, p) # Compute the Hopf point
    t = 0.0
    u = zeros(dimension(eq))
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    s = nlsolve(p0 -> maximum(real(eigen(jeq(u, p0, t)).values)) - 1e-11, p)
    return (mu=s.zero[1],eig=eigen(jeq(u, s.zero[1], t)).values)
end

function set_params_1()
    return [0.15, -0.5, 1.204, 0.24, 0.5628, 14.5756, 0.1724, 54.1162, 3.5294e+03, 5.3, 16.9]
end          #b,     a,     ρ,   x_α,   c_α,     c_h,    I_α,    k_α,      k_h,      m, m_T


function set_params_2() 
    return [0.15, -0.5, 1.204, 0.24, 0.9426, 14.5756, 0.1724, 61.3039, 3.3183e+03, 5.3, 16.9]
end         #b,     a,     ρ,   x_α,   c_α,     c_h,    I_α,    k_α,      k_h,      m, m_T
   
function flutter_eq_RSPA_model1(u, p, t) #  Equation of motion of Flutter model with CBC
    U = p[1] 
    Model=RSPA_model(Continuation.set_params_1(),U)
    J=Model.J
    M2=Model.M2
    du = J * u
       
    ka2=750.67;ka3=5007.38
    nonlinear_h = -M2[1, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)
    nonlinear_theta = -M2[2, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)

    du[2] += nonlinear_h 
    du[4] += nonlinear_theta 
    return du
end

function flutter_eq_RSPA_model2(u, p, t) #  Equation of motion of Flutter model with CBC_RSPA model-2
    U = p[1] 
    Model=RSPA_model(Continuation.set_params_2(),U)
    J=Model.J
    M2=Model.M2
    du = J * u
       
    #ka2=829.50;ka3=4037.19  
    ka2=711.88;ka3=3222.52    
    nonlinear_h = -M2[1, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)
    nonlinear_theta = -M2[2, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)

    du[2] += nonlinear_h 
    du[4] += nonlinear_theta 
    return du
end

function flutter_eq_RSPA_model1_u(u, p, t) #  Equation of motion of Flutter model with CBC
    U = p[1] 
    Model=RSPA_model(Continuation.set_params_1(),U)
    J=Model.J
    M2=Model.M2
    du = J * u
       
    ka2=751.18; ka3=4769.10
    nonlinear_h = -M2[1, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)
    nonlinear_theta = -M2[2, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)

    du[2] += nonlinear_h 
    du[4] += nonlinear_theta 
    return du
end

function RSPA_model(p,U) #  Equation of motion of Flutter model with CBC
    b,a,ρ,x_α,c_α,c_h,I_α,k_α,k_h,m,m_T=p
    
    K1=[-k_h*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2))) -(k_α - 0.9999999999999999π*ρ*(U^2)*(0.5 + a)*(b^2))*((-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2))) - 0.9999999999999999b*π*ρ*(U^2)*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2))) 0.013649999999999999b*π*ρ*(U^3)*(0.5 + a)*((-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2))) - 0.013649999999999999π*ρ*(U^3)*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2))); -k_h*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) -(k_α - 0.9999999999999999π*ρ*(U^2)*(0.5 + a)*(b^2))*(true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) - 0.9999999999999999b*π*ρ*(U^2)*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) 0.013649999999999999b*π*ρ*(U^3)*(0.5 + a)*(true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) - 0.013649999999999999π*ρ*(U^3)*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))); -0.0 U / b (-0.013649999999999999(U^2)) / (b^2)]

    D1=[0.9999999999999999U*π*ρ*(0.5 + a)*(b^2)*((-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2))) - (c_h + 0.9999999999999999U*b*π*ρ)*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T 
+ π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2))) -(c_α + U*π*ρ*(b^3)*(0.5 - 0.9999999999999999a)*(0.5 - a))*((-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2))) - U*π*ρ*(b^2)*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2)))*(1.5 - 0.9999999999999999a) 0.216015π*ρ*(U^2)*(0.5 + a)*(b^2)*((-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2))) - 0.216015b*π*ρ*(U^2)*((true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2))); 0.9999999999999999U*π*ρ*(0.5 + a)*(b^2)*(true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) - (c_h + 0.9999999999999999U*b*π*ρ)*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) -(c_α + U*π*ρ*(b^3)*(0.5 - 0.9999999999999999a)*(0.5 - a))*(true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) - U*π*ρ*(b^2)*(1.5 - 0.9999999999999999a)*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) 0.216015π*ρ*(U^2)*(0.5 + a)*(b^2)*(true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) - 0.216015b*π*ρ*(U^2)*((-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))); 1.0 / b 0.5 - a (-0.3455U) / b]

    J1 = [0 1 0 0 0 0]
    J2 = [K1[1, 1] D1[1, 1] K1[1, 2] D1[1, 2] K1[1, 3] D1[1, 3]]
    J3 = [0 0 0 1 0 0]
    J4 = [K1[2, 1] D1[2, 1] K1[2, 2] D1[2, 2] K1[2, 3] D1[2, 3]]
    J5 = [0 0 0 0 0 1]
    J6 = [K1[3, 1] D1[3, 1] K1[3, 2] D1[3, 2] K1[3, 3] D1[3, 3]]

    J = [J1; J2; J3; J4; J5; J6]

    M2=[(true + (((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))*(b*m*x_α - a*π*ρ*(b^3))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2))) / (m_T + π*ρ*(b^2)) (-((b*m*x_α - a*π*ρ*(b^3)) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)))) / (m_T + π*ρ*(b^2)) 0; (-((b*m*x_α - a*π*ρ*(b^3)) / (m_T + π*ρ*(b^2)))) / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)) true / (I_α + (-((b*m*x_α - a*π*ρ*(b^3))^2)) / (m_T + 
    π*ρ*(b^2)) + π*ρ*(b^4)*(0.125 + a^2)) 0; 0.0 0.0 1.0]
    
    return (J=J,M2=M2)
end


dimension(::typeof(flutter_eq_RSPA_model1_u)) = 6
dimension(::typeof(flutter_eq_RSPA_model1)) = 6
dimension(::typeof(flutter_eq_RSPA_model2)) = 6

function Flutter_solve(n, p, ka2, ka3, p_i, eq, p_r2) 
    # Flutter_solve is function that solves periodic boundary value problem of flutter with using normal form as an initial guess
    hp=Hopf_point(eq, [p_i])
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    J=jeq(zeros(6),[hp],0.0)
    eig_v=eigvals(J)
    re_eig=abs.(real.(eig_v))
    ω=abs(imag(eig_v[argmin(re_eig)]))
    ini_T = 2π / ω
    t = range(0, ini_T, length=n + 1)[1:end-1]
    eig=eigvecs(J)[:,argmin(re_eig)]
    arg = atan(imag(eig[1]) / real(eig[1]))
    eig = exp((-arg + π / 2) * 1im) * eig

    R2 = p_r2[1] * (p[1] - hp) / (p_r2[2]*ka3 + p_r2[3]*ka2^2)
    R = sqrt(R2)
    u1 = 2*real([R * eig[1] * exp.(1im * ω * t)'; R * eig[2] * exp.(1im * ω * t)'; R * eig[3] * exp.(1im * ω * t)'; R * eig[4] * exp.(1im * ω * t)'; R * eig[5] * exp.(1im * ω * t)'; R * eig[6] * exp.(1im * ω * t)'])
    u1=vec(u1)

    z_tol = 1e-9
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    ee=1e-4
    s=LCO_solve(u1,ini_T,eq,jeq,p,ee,z_tol)
    return s
end

function continuation_flutter(N, sp,  ds, eq, ka2, ka3,p_i,p_r2) #Get stable LCO from the ODE numerical solutions- initial collocation points
    ee = 1e-4
    tol = 1e-7
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    jeq2 = (u, p, t) -> ForwardDiff.jacobian(p -> eq(u, p, t), p)[:, 1]
    
    hp=Hopf_point(eq, [p_i])
   
    s=Flutter_solve(N, [hp-0.0001], ka2, ka3, p_i, eq, p_r2)
    s2=Flutter_solve(N, [hp-0.0002], ka2, ka3, p_i, eq, p_r2)

    du = s2.u - s.u
    du = [du; -0.0001]
    dv = du / norm(du)
    dim = 6
    ll = dim * N + 1
    V = [zeros(ll) for _ in 1:sp]
    P = zeros(sp)
    mu = zeros(Complex,dimension(eq)*N,sp)
    mu_m = zeros(sp)
    pu = s2.u
    pu = [pu; hp-0.0002]
    V[1] = vec(s2.u)
    P[1] = hp-0.0002

    uu = pu + ds * dv
    μ = uu[end]
    p = [μ]
    u = uu[1:end-1]
    u1 = u[1:end-1]
    ini_T = u[end]
    z_tol = tol * 0.1
    M=[0]
#    M=Monodromy_compute(eq, s2.u, [P[1]], N, [1,2])
    mu_m[1]=maximum(abs.(M))
  #  mu[:,1]=M
    for i in 2:sp
        s2 = LCO_solve2(u1, ini_T, eq, jeq, jeq2, p, ee, z_tol, dv, pu, ds)
        s2.err
        du = s2.u - pu
        dv = du / norm(du)
        V[i] = s2.u[1:end-1]
        P[i] = s2.u[end]
        pu = s2.u
        uu = pu + ds * dv
        p = [uu[end]]
        u = uu[1:end-1]
        u1 = u[1:end-1]
        ini_T = u[end]
        z_tol = tol * 0.1
#        M=Monodromy_compute(eq, s2.u[1:end-1], [s2.u[end]], N, [1,2])
        mu_m[i]=maximum(abs.(M))
 #       mu[:,i]=M
    end
    return (V=V, P=P,mu=mu,mu_m=mu_m)
end

function flutter_eq_RSPA_model_update(u, p, t) #  Equation of motion of Flutter model with CBC
    U = p[1] 
    Model=RSPA_model(Continuation.set_params_1(),U)
    J=Model.J
    M2=Model.M2
    du = J * u
       
    ka2=p[2];ka3=p[3]
    nonlinear_h = -M2[1, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)
    nonlinear_theta = -M2[2, 2] * (ka2 * u[3]^2 + ka3 * u[3]^3)

    du[2] += nonlinear_h 
    du[4] += nonlinear_theta 
    return du
end

dimension(::typeof(flutter_eq_RSPA_model_update)) = 6

function model_update_solve(p, uu, eq)   
    z_tol = 1e-9
    jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
    ee=1e-4
    u1=uu[1:end-1]
    ini_T=uu[end]
    s=LCO_solve(u1,ini_T,eq,jeq,p,ee,z_tol)
    return s
end

end