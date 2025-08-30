## ___bb-20250720-snapshot.py
## % time python ___bb-20250720-snapshot.py seed e_L e_Dr e_v0 e_Dth &

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def theta_in2pi(theta):
	while theta<0.0:
		theta+=2.0*np.pi
	while theta>=2.0*np.pi:
		theta-=2.0*np.pi
	return theta

seed=int(sys.argv[1])
e_L=int(sys.argv[2])
e_Dr=int(sys.argv[3])
e_v0=float(sys.argv[4])
e_Dth=int(sys.argv[5])

if e_v0!=3.5:
	data_parameter=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
elif e_v0==3.5:
	data_parameter=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
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

tau_snapshot=Tau
#tau_snapshot=500000

if e_v0!=3.5:
	data_movie=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=2))
elif e_v0==3.5:
	data_movie=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=2))

fig=plt.figure(figsize=(10,8))
ax=[fig.add_subplot(1,4,(1,3))]
fs=10
lw=1
ms=4/L_square
ds=3#*L_square
cmap=plt.get_cmap("hsv")
ax[0].set_xticks([0,L_square])
ax[0].set_xticklabels(["0.0","%.1f"%L_square],fontsize=fs)
ax[0].set_xlim(0,L_square)
ax[0].set_yticks([0,L_square])
ax[0].set_yticklabels(["0.0","%.1f"%L_square],fontsize=fs)
ax[0].set_ylim(0,L_square)
ax[0].set_aspect(1)

j=0
J_snapshot=0
while j<len(data_movie) and J_snapshot==0:
	tau=int(data_movie[j,0])
	M_now=int(data_movie[j,1])
	g_max=float(data_movie[j,2])
	S1=float(data_movie[j,3])
	j+=1
	if tau==tau_snapshot:
		ax[0].set_title("$t=$%.2f, $M=$%d, $g_{max}=$%.2f, $S=$%.2f"%(tau*dt,M_now,g_max,S1),fontsize=fs)
		for i in range(M_now):
			x_i=float(data_movie[j+i,1])
			y_i=float(data_movie[j+i,2])
			th_i=float(data_movie[j+i,3])
			ax[0].plot(x_i,y_i,marker="o",markersize=ms,color=cmap(theta_in2pi(th_i)*0.5/np.pi))
			ax[0].plot([x_i,x_i+0.01*ds*np.cos(th_i)],[y_i,y_i+0.01*ds*np.sin(th_i)],linewidth=lw,color=cmap(theta_in2pi(th_i)*0.5/np.pi))
		J_snapshot=1
	j+=M_now

imageForColorbar=ax[0].imshow(np.array([[]]),cmap="hsv",vmin=0.0,vmax=2.0*np.pi,alpha=1.0)
cbar=plt.colorbar(imageForColorbar,ax=ax[0],shrink=0.8,location="right",orientation="vertical")
cbar.set_ticks(np.linspace(0,2*np.pi,5))
cbar.set_ticklabels(["0",r"0.5$\pi$",r"$\pi$",r"1.5$\pi$",r"2$\pi$"])
cbar.ax.tick_params(labelsize=fs)

ax[0].plot([],[],label="seed$=$%d"%seed)
ax[0].plot([],[],label="$D_r=$%.3e"%Dr)
ax[0].plot([],[],label="$v_0=$%.3e"%v0)
ax[0].plot([],[],label=r"$D_\theta=$%.0e"%Dth)
ax[0].plot([],[],label="$dt=$%.0e"%dt)
ax[0].plot([],[],label="$T=$%d"%Tau)
ax[0].plot([],[],label="$T_{relax}=$%d"%Tau_relax)
ax[0].plot([],[],label="interval$=$%d"%interval)
ax[0].plot([],[],label="$L=$%.1f"%L_square)
ax[0].plot([],[],label="$R=$%.1f"%R_inter)
ax[0].plot([],[],label="$N_s=$%d"%Ns)
ax[0].plot([],[],label="$p_0=$%.2f"%p0)
ax[0].plot([],[],label="$q_0=$%.2f"%q0)
ax[0].plot([],[],label="$M_0=$%d"%M0)
ax[0].plot([],[],label="$dr=$%.0e"%dr)
ax[0].plot([],[],label="$N_r=$%d"%Nr)
ax[0].legend(loc="upper left",bbox_to_anchor=(1.30,0.00,1,1),ncol=1,handlelength=0,handletextpad=0,fontsize=fs)

if e_v0!=3.5:
	fig.savefig("__bb-20250720-snapshot%d-%d-%d-%d-%d-%d.pdf"%(tau_snapshot,seed,e_L,e_Dr,e_v0,e_Dth))
elif e_v0==3.5:
	fig.savefig("__bb-20250720-snapshot%d-%d-%d-%d-%.1f-%d.pdf"%(tau_snapshot,seed,e_L,e_Dr,e_v0,e_Dth))

