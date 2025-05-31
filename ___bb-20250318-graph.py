## ___bb-20250318-graph.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def my_round(c):
    return int(c)+(c-int(c)>=0.5)-(c-int(c)<-0.5)

seed=int(sys.argv[1])
e_L=int(sys.argv[2])
e_Dr=int(sys.argv[3])
e_v0=int(sys.argv[4])
e_Dth=int(sys.argv[5])

data_parameter=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
data_movie=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=1))
data_diagram=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1+len(data_movie)))

j=0
while j<len(data_movie):
    if j==0:
        data_graph=data_movie[:1]
    else:
        data_graph=np.append(data_graph,data_movie[j:j+1],axis=0)
    j+=1+int(data_movie[j,1])

#data_parameter=np.array(pd.read_csv("_bb-20250318-data_parameter-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python"))
Dr=float(data_parameter[0,1])
v0=float(data_parameter[0,2])
Dth=float(data_parameter[0,3])
dt=float(data_parameter[0,4])
Tau=int(data_parameter[0,5])
Tau_relax=int(data_parameter[0,6])
interval=int(data_parameter[0,7])
L_square=float(data_parameter[0,8])
R_inter=float(data_parameter[0,9])
Ns=int(data_parameter[0,10])
p0=float(data_parameter[0,11])
q0=float(data_parameter[0,12])
M0=int(data_parameter[0,13])
dr=float(data_parameter[0,14])
Nr=int(data_parameter[0,15])

#data_graph=np.array(pd.read_csv("_bb-20250318-data_graph-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python"))
list_M=np.zeros(shape=my_round(Tau/interval)+1)
list_rho=np.zeros(shape=my_round(Tau/interval)+1)
list_gmax=np.zeros(shape=my_round(Tau/interval)+1)
list_S=np.zeros(shape=my_round(Tau/interval)+1)
list_s=np.zeros(shape=my_round(Tau/interval)+1)
list_Nf=np.zeros(shape=my_round(Tau/interval)+1)
for j in range(len(data_graph)):
    list_M[j]=float(data_graph[j,1])
    list_rho[j]=float(data_graph[j,1])/(L_square**2)
    list_gmax[j]=float(data_graph[j,2])
    list_S[j]=float(data_graph[j,3])
    list_s[j]=float(data_graph[j,4])
    list_Nf[j]=float(data_graph[j,5])

#data_diagram=np.array(pd.read_csv("_bb-20250318-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python"))
#M_=float(data_diagram[0,0])
#rho_=float(data_diagram[0,0])/(L_square**2)
S_=float(data_diagram[0,1])
#s_=float(data_diagram[0,2])
list_g=np.zeros(shape=Nr)
g_r2=1.0
r2=0.0
J_r2=0
for j in range(Nr):
    list_g[j]=float(data_diagram[0,3+j])
    if J_r2==0 and list_g[j]<1.0:
        J_r2=1
    elif J_r2==1 and g_r2<=list_g[j]:
        g_r2=list_g[j]
        r2=j*dr
g_max=float(data_diagram[0,-1])

ax_v=10

fig=plt.figure(figsize=(12,8))
ax=[fig.add_subplot(2,ax_v,(1,ax_v-1)),fig.add_subplot(2,ax_v,(ax_v+1,2*ax_v-1))]
ax.append(ax[0].twinx())
fs=8
lw=6
tl=10
tw=6

ax[0].set_xlabel("$r$",labelpad=0,fontsize=fs)
ax[0].set_xticks(np.arange(0.0,0.5*L_square+R_inter,R_inter))
ax0_xticklabels=["0","$R$"]
ax0_xticklabels.extend(["%d$R$"%my_round(j/R_inter) for j in np.arange(R_inter*2,L_square/2,R_inter)])
ax0_xticklabels.append("$L/2$")
ax[0].set_xticklabels(ax0_xticklabels,fontsize=fs)
ax[0].set_xlim(0,R_inter*5)
ax[0].set_ylabel("$g$",labelpad=0,fontsize=fs)
ax0_yticks=np.arange(0,2.1,0.5)
ax[0].set_yticks(ax0_yticks)
ax[0].set_yticklabels(["%.1f"%j for j in ax0_yticks],fontsize=fs)
ax[0].set_ylim(0,2)
ax[0].tick_params(axis="both",right=False,length=tl,width=tw)
#ax[0].axhline(y=1.0,linestyle="--",color="k",linewidth=lw)
ax[0].hlines(y=g_max,xmin=0.0,xmax=r2,linestyle="--",color="k",linewidth=lw)
ax[0].plot(np.arange(0.0,0.5*L_square,dr),list_g,color="g",linewidth=lw)
#ax[0].set_title("$g_{max}=$%.2f"%g_max,fontsize=fs)

ax[1].set_xlabel("$t$",labelpad=0,fontsize=fs)
ax[1].set_xticks(np.linspace(0,Tau,5)*dt)
ax[1].set_xticklabels(["%.1e"%j for j in np.linspace(0,Tau,5)*dt],fontsize=fs)
ax[1].set_xlim(0.0,Tau*dt)
ax[1].set_ylabel("$S$",labelpad=0,fontsize=fs)
ax[1].set_yticks(np.linspace(0,1,5))
ax[1].set_yticklabels(["%.2f"%j for j in np.linspace(0,1,5)],fontsize=fs)
ax[1].set_ylim(0,1)
ax[1].tick_params(axis="both",right=True,length=tl,width=tw)
ax[1].plot(np.arange(0,Tau+1,interval)*dt,list_S,color="b",linewidth=lw)

ax[2].tick_params(axis="y",right=False,labelright=False)
ax[2].plot([],[],label="seed$=$%d"%seed)
ax[2].plot([],[],label="$D_r=$%.0e"%Dr)
ax[2].plot([],[],label="$v_0=$%.0e"%v0)
ax[2].plot([],[],label=r"$D_\theta=$%.0e"%Dth)
ax[2].plot([],[],label="$dt=$%.0e"%dt)
ax[2].plot([],[],label="$T=$%d"%Tau)
ax[2].plot([],[],label="$T_{relax}=$%d"%Tau_relax)
ax[2].plot([],[],label="interval$=$%d"%interval)
ax[2].plot([],[],label="$L=$%.1f"%L_square)
ax[2].plot([],[],label="$R=$%.1f"%R_inter)
ax[2].plot([],[],label="$N_s=$%d"%Ns)
ax[2].plot([],[],label="$p_0=$%.2f"%p0)
ax[2].plot([],[],label="$q_0=$%.2f"%q0)
ax[2].plot([],[],label="$M_0=$%d"%M0)
ax[2].plot([],[],label="$dr=$%.0e"%dr)
ax[2].plot([],[],label="$N_r=$%d"%Nr)
ax[2].legend(loc="upper left",bbox_to_anchor=(1.05,-0.05,1,1),ncol=1,handlelength=0,handletextpad=0,fontsize=fs)

fig.savefig("__bb-20250318-graph_draft-%d-%d-%d-%d-%d.png"%(seed,e_L,e_Dr,e_v0,e_Dth))
