module Flutter_opt

using StaticArrays, NLsolve, NLopt

function measure_error(ka2::Real, ka3::Real) # Computing the numerical prediction error (Model1)
    # measure_error function gives the error=[numerical_amplitude-measured_amplitude] at 2 measured wind speeds
    us = Vector{Float64}(undef, 3)
    mh = Vector{Float64}(undef, 3)
    # Measuremen results
    us[1] = 17.2
    us[2] = 16.5
    us[3] = 15.7

    mh[1] = 0.007523
    mh[2] = 0.010891
    mh[3] = 0.0117

    us=us-ones(3)*17.9286 # Normalise control parameter

    e=0

    v2=us
    R2=-0.0087943864526163*v2/(-4.55040566767517e-05*ka3 + 8.41450867929661e-07*ka2^2)
    R=sqrt.(R2)
    eig = -0.0155041+0.0159816im
    h_portion=abs(eig)*2
    hn=R*h_portion
    for i=1:length(hn)
        e+=(hn[i]-mh[i])^2
    end
    return e*1000 # Scale of measured error is mm
end

function measure_error2(ka2::Real, ka3::Real) # Computing the numerical prediction error (Model2)
    # measure_error function gives the error=[numerical_amplitude-measured_amplitude] at 2 measured wind speeds
    us = Vector{Float64}(undef, 4)
    mh = Vector{Float64}(undef, 4)
    # Measuremen results
    us[1] = 25.8
    us[2] = 25.0
    us[3] = 24.1
    us[4] = 23.3

    us=us-ones(4)*26.1745 # Normalise control parameter

    mh[1] = 0.0088
    mh[2] = 0.0151
    mh[3] = 0.0224
    mh[4] = 0.0258
    e=0

    v2=us

    R2=-0.00726461306670313*v2/(-3.21407245501212e-05*ka3 + 4.97046686878214e-07*ka2^2)
    R=sqrt.(R2)
    eig=0.022824369459937155 - 0.025183235928075854im
    h_portion=abs(eig)*2
    hn=R*h_portion

    for i=1:length(hn)
        e+=(hn[i]-mh[i])^2
    end
    return e*1000  # Scale of measured error is mm
end

function my_fun(x::Vector, grad::Vector) # Function to minimize (Model 1)
    a = x[1]
    b = x[2]

    if length(grad) > 0
        h = 0.0001
        grad[1] = (measure_error(a + h, b) - measure_error(a, b)) / h
        grad[2] = (measure_error(a, b + h) - measure_error(a, b)) / h
    end
    d = (measure_error(a, b))
    return d
end

function my_fun2(x::Vector, grad::Vector) # Function to minimize (Model 2)
    a = x[1]
    b = x[2]

    if length(grad) > 0
        h = 0.0001
        grad[1] = (measure_error2(a + h, b) - measure_error2(a, b)) / h
        grad[2] = (measure_error2(a, b + h) - measure_error2(a, b)) / h
    end
    d = (measure_error2(a, b))
    return d
end


function my_con(x::Vector, grad::Vector) # Constraint of the optimization (Model 1)
    if length(grad) > 0
        grad[1] = -2 * 8.41659973544257e-07 * x[1]
        grad[2] = 4.52216831510469e-05
    end
    -((-4.55040566767517e-05 * (x[2])) + 8.41450867929661e-07 * (x[1])^2)
end

function my_con2(x::Vector, grad::Vector) # Constraint of the optimization (Model 2)
    if length(grad) > 0
        grad[1] = -2 * 4.97046686878214e-07 * x[1] 
        grad[2] = 3.21407245501212e-05
    end
    -(-3.21407245501212e-05*x[2] + 4.97046686878214e-07*x[1]^2)
end

#Optimization using NLopt package
function optimize_NLstiff(x1, x2) # Model 1
    opt = Opt(:LD_MMA, 2)
    opt.lower_bounds = [100, 100]
    opt.upper_bounds = [1000, 7000]
    opt.xtol_abs = 1.5 * 1e-2
    opt.min_objective = my_fun
    inequality_constraint!(opt, (x, g) -> my_con(x, g), 1e-3)
    opt = optimize(opt, [x1, x2])
    return opt
end

function optimize_NLstiff2(x1, x2) # Model 2
    opt = Opt(:LD_MMA, 2)
    opt.lower_bounds = [100, 100]
    opt.upper_bounds = [1000, 7000]
    opt.xtol_abs = 1e-2
    opt.min_objective = my_fun2
    inequality_constraint!(opt, (x, g) -> my_con2(x, g), 1e-3)
    opt = optimize(opt, [x1, x2])
    return opt
end

end #modlue
