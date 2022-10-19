using HopfNormalForm
using Plots
using PGFPlotsX
using LaTeXStrings



sys1=Systems.example1()
sys2=Systems.example2()

sys1.LCO_R2
sys2.LCO_R2

p_r2_1=[-0.0087943864526163,-4.55040566767517e-05,8.41450867929661e-07]
p_r2_2=[-0.00726461306670313,-3.21407245501212e-05,4.97046686878214e-07]

sys1.SLC
sys2.SLC

## Model 1
m1=Flutter_opt.optimize_NLstiff(751.0,5006.7)
ka2=750.67;ka3=5007.38
eq=Continuation.flutter_eq_RSPA_model1
p_i=18.0

#Continuation.Hopf_point2(eq,[p_i]).eig
flutter_freq=15.4789/2/pi


N=15
cont=Continuation.continuation_flutter(N, 250,0.04,eq, ka2, ka3,p_i,p_r2_1)
amp=Continuation.amp_LCO(cont.V, 1, N)
pp=cont.P
pp=pp-ones(length(pp))*17.9572

Vel3=[14.9, 15.6, 16.5, 17.3, 14.9, 15.7, 16.5, 17.2]
h3=[ 0.0196, 0.0251, 0.0322, 0.0448, 0.0158, 0.0117, 0.0109, 0.0078]
Vel3=Vel3-ones(length(Vel3))*17.9286 # Normalise the measured wind speed
#scatter!(Vel3,h3)
#plot(pp,amp/2)
#scatter!(Vel3,h3)

bd = @pgf Axis(
    {
        xlabel = L"$U_f-U$",
        ylabel = "Heave amplitude (m)",
        legend_pos = "north west",
        height = "7cm",
        width = "11cm",
        ymin = 0,
        ymax = 9e-2,
        mark_options = {scale = 1.5},
    },
    Plot({color = "blue"}, Coordinates(pp, amp)),
    Plot({color = "red", only_marks}, Coordinates(Vel3, h3*2)),
)
@pgf bd["every axis title/.style"] = "below right,at={(0,1)}";
bd
##
pgfsave("./Figure/bd_m1.pdf", bd)

## Model 2

m2=Flutter_opt.optimize_NLstiff2(720,3200.7)

ka2=711.88;ka3=3222.52
eq=Continuation.flutter_eq_RSPA_model2
#Continuation.Hopf_point2(eq,[p_i]).eig
#flutter_freq=15.26/2/pi
p_i=23.0
#cont=Continuation.continuation_flutter(30, 500,0.04,eq, ka2, ka3,p_i,p_r2_2)

cont=Continuation.continuation_flutter(30, 400,0.05,eq, ka2, ka3,p_i,p_r2_2)
amp=Continuation.amp_LCO(cont.V, 1, 30)
pp=cont.P
pp=pp-ones(length(pp))*23.6951 # Normalise the control parameter

Vel3=[20.7, 21.5, 22.4, 23.2, 20.6, 21.6, 22.4, 23.3, 24.1, 25.0, 25.8]
Vel3=Vel3-ones(length(Vel3))*26.1745 # Normalise the measured wind speed
h3=[ 0.0523, 0.0645,0.0754,0.0875,0.0341,0.0344,0.0288,0.0258, 0.0224, 0.0151, 0.0088 ]

plot(pp,amp/2)
scatter!(Vel3,h3)

bd2 = @pgf Axis(
    {
        xlabel = L"$U_f-U$",
        ylabel = "Heave amplitude (m)",
        legend_pos = "north west",
        height = "7cm",
        width = "11cm",
        ylabel_shift="-5pt",
        ymin = 0,
        ymax = 20e-2,
        ytick=[0.1, 0.2],
        mark_options = {scale = 1.5},
    },
    Plot({color = "blue"}, Coordinates(pp, amp)),
    Plot({color = "red", only_marks}, Coordinates(Vel3, h3*2)),
)
@pgf bd2["every y tick scale label/.style"] = "at={(yticklabel cs:1)},anchor=south west";

bd2
pgfsave("./Figure/bd_m2.pdf", bd2)



