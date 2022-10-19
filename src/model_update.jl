using HopfNormalForm
using Plots
using PGFPlotsX
using LaTeXStrings
using LinearAlgebra

p_r2_1=[-0.0087943864526163,-4.55040566767517e-05,8.41450867929661e-07]
p_r2_2=[-0.00726461306670313,-3.21407245501212e-05,4.97046686878214e-07]

ka2=750.67;ka3=5007.38
eq=Continuation.flutter_eq_RSPA_model1
p_i=18.0
N=15
cont=Continuation.continuation_flutter(N, 250,0.04,eq, ka2, ka3,p_i,p_r2_1)
amp=Continuation.amp_LCO(cont.V, 1, N)
pp=cont.P

# Find solutions at the measured wind speed using collocation

sign1=zeros(length(pp))

for ii=2:length(pp)
    sign1[ii]=pp[ii]-pp[ii-1]
end
# sign1[112] Fold point index 111

Vel3=[14.9, 15.6, 16.5, 17.3, 14.9, 15.7, 16.5, 17.2]
Vel3_s=[14.9, 15.6, 16.5, 17.3]
Vel3_u=[14.9, 15.7, 16.5, 17.2]
s_i=zeros(4)
u_i=zeros(4)

for ii=1:length(Vel3_s)
    zz=(pp[112:end]-Vel3_s[ii]*ones(length(pp[112:end])))
    zz=broadcast(abs, zz)
    s_i[ii]=argmin(zz)+111
end

for ii=1:length(Vel3_s)
    zz=(pp[1:111]-Vel3_u[ii]*ones(length(pp[1:111])))
    zz=broadcast(abs, zz)
    u_i[ii]=argmin(zz)
end

uu1=[cont.V[Int(s_i[ii])] for ii=1:4]
uu2=[cont.V[Int(u_i[ii])] for ii=1:4]

function cost_fun_update(ka2,ka3,uu1,uu2) # Define a cost function for optimisation
amp_=zeros(4)
amp_2=zeros(4)

    for ii=1:length(s_i)
        p=[Vel3_s[ii],ka2,ka3]
        uu=uu1[ii]
        eq_=Continuation.flutter_eq_RSPA_model_update
        s=Continuation.model_update_solve(p, uu, eq_)  
        g = Continuation.get_sol(s.u, 6, N, 1, 3)
        uu = g.u[:, 1]
        uu1[ii]=s.u
        amp_[ii] = maximum(uu) - minimum(uu)
    end

    for ii=1:length(u_i)
        p=[Vel3_u[ii],ka2,ka3]
        uu=uu2[ii]
        eq_=Continuation.flutter_eq_RSPA_model_update
        s=Continuation.model_update_solve(p, uu, eq_)  
        g = Continuation.get_sol(s.u, 6, N, 1, 3)
        uu = g.u[:, 1]
        uu2[ii]=s.u
        amp_2[ii] = maximum(uu) - minimum(uu)
    end

return    (c=norm(vcat(amp_,amp_2)-h3*2),uu1=uu1,uu2=uu2)
end

function mod_up_J(ka2,ka3,uu1,uu2)
    epsil=1e-1
    cc=cost_fun_update(ka2,ka3,uu1,uu2)
    cc1=cost_fun_update(ka2+epsil,ka3,uu1,uu2)
    cc2=cost_fun_update(ka2,ka3+epsil,uu1,uu2)
    d1=(cc1.c-cc.c)/epsil
    d2=(cc2.c-cc.c)/epsil
return    (d1=d1,d2=d2)
end

function model_update(ka2,ka3,uu1,uu2,iter,delta)
    c2=100
    c1=0
    k=0
    for ii=1:iter
        dir=mod_up_J(ka2,ka3,uu1,uu2)
        ka2=ka2-dir.d1*delta
        ka3=ka3-dir.d2*delta
        cc=cost_fun_update(ka2,ka3,uu1,uu2)
        uu1=cc.uu1
        uu2=cc.uu2
        c1=cc.c
        k=ii
        if c1>c2
            break
        else
            c2=c1
        end    
    end
    return (c=c1,ka2=ka2,ka3=ka3,k=k,uu1=uu1,uu2=uu2)
end

iter=1000
ka2=751.18; ka3=4769.10
uu1=[cont.V[Int(s_i[ii])] for ii=1:4]
uu2=[cont.V[Int(u_i[ii])] for ii=1:4]
delta=500
mu=model_update(ka2,ka3,uu1,uu2,iter,delta)

c0=cost_fun_update(mu.ka2,mu.ka3,mu.uu1,mu.uu2)

iter=10000
delta=1000
mu=model_update(mu.ka2,mu.ka3,mu.uu1,mu.uu2,iter,delta)
norm([mod_up_J(mu.ka2,mu.ka3,mu.uu1,mu.uu2).d1,mod_up_J(mu.ka2,mu.ka3,mu.uu1,mu.uu2).d2])


ka2=751.18; ka3=4769.10
eq=Continuation.flutter_eq_RSPA_model1_u
p_i=18.0
N=15
cont_u=Continuation.continuation_flutter(N, 280,0.04,eq, ka2, ka3,p_i,p_r2_1)
amp_u=Continuation.amp_LCO(cont_u.V, 1, N)
pp_u=cont_u.P



plot(pp,amp)
plot!(pp_u,amp_u)
scatter!(Vel3,h3*2)
h3=[ 0.0196, 0.0251, 0.0322, 0.0448, 0.0158, 0.0117, 0.0109, 0.0078]

bd = @pgf Axis(
    {
        xlabel = "Wind speed (m/s)",
        ylabel = "Heave amplitude (m)",
        legend_pos = "north west",
        height = "7cm",
        width = "11cm",
        ymin = 0,
        ymax = 9e-2,
        xmax=19.0,
        mark_options = {scale = 1.5},
    },
    Plot({color = "blue"}, Coordinates(pp, amp)),
    Plot({color = "red"}, Coordinates(pp_u, amp_u)),
    Plot({color = "red", only_marks}, Coordinates(Vel3, h3*2)),
)
@pgf bd["every axis title/.style"] = "below right,at={(0,1)}";
bd

pgfsave("./Figure/bd_m1_upd.pdf", bd)




bd2 = @pgf Axis(
    {
        xlabel = "Wind speed (m/s)",
        ylabel = "Heave amplitude (m)",
        legend_pos = "north west",
        height = "7cm",
        width = "11cm",
        ymin = 0,
        ymax = 9e-2,
        xmax=19.0,
        mark_options = {scale = 1.5},
    },
    Plot({color = "blue"}, Coordinates(pp, amp)),
    Plot({color = "red", only_marks}, Coordinates(Vel3, h3*2)),
)
bd2

pgfsave("./Figure/bd_m1_2.pdf", bd2)