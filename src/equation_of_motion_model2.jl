using SymEngine
using LinearAlgebra
using ZChop
using ForwardDiff
using HopfNormalForm

function swapcol!(x,i,j)
    for k in axes(x,1) # edited according to next answer
      idata = x[k,i]
      x[k,i] = x[k,j]
      x[k,j] = idata
    end
end

Uf=Continuation.Hopf_point(Continuation.flutter_eq_RSPA_model2, [22.0]) #Compute the Hopf point

U=Uf
eq=Continuation.flutter_eq_RSPA_model2
jeq = (u, p, t) -> ForwardDiff.jacobian(u -> eq(u, p, t), u)
J=jeq(zeros(6),[Uf],0.0)

e_J=eigen(J) # Compute eigenvalue at the Hopf point
eiVec=e_J.vectors

##Re arrange eigenvectors in the form of HopfNormalForm.jl
swapcol!(eiVec,6,1)
swapcol!(eiVec,5,2)
swapcol!(eiVec,5,6)
swapcol!(eiVec,3,4)
##
eJ=inv(eiVec)*J*eiVec
eJ=zchop(eJ,1e-9) # Chop out numerical garbage

v=[symbols("v_$i") for i in 1:7]
@vars ka2, ka3

U=Uf+v[1] # set U with folding parameter

b,a,ρ,x_α,c_α,c_h,I_α,k_α,k_h,m,m_T=Continuation.set_params_2()
model=Continuation.RSPA_model(Continuation.set_params_2(),U)
J=model.J
M2=model.M2
pJ=inv(eiVec)*J*eiVec

for i=1:6
    for j=1:6
    pJ[i,j]=expand(HopfNormalForm.chop_expr(pJ[i,j],1e-9))
    end
end

original_coordinate=eiVec*v[2:7]
u_3=original_coordinate[3]
u_3=expand(HopfNormalForm.chop_expr(u_3,1e-9))

nonlinear_h = -M2[1, 2] * (ka2 * u_3^2 + ka3 * u_3^3)
nonlinear_theta = -M2[2, 2] * (ka2 * u_3^2 + ka3 * u_3^3)

non_stiff=[0,nonlinear_h,0,nonlinear_theta,0,0]
non_stiff=inv(eiVec)*non_stiff

for j=1:6
    non_stiff[j]=HopfNormalForm.chop_expr(expand(non_stiff[j]),1e-9)
end

EM=pJ*v[2:7]+non_stiff
for i=1:6
    EM[i]=HopfNormalForm.chop_expr(expand(EM[i]),1e-9)
end

## Print below and change v_i->v[i] for i=1:7. 
EM[6]