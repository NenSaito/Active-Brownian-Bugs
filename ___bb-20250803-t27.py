## ___bb-20250803-t27.py
## % time python ___bb-20250803-t27.py seed e_L n_Dr n_v0 e_Dth &

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

def my_round(c):
	return int(c)+int(c-int(c)>=0.5)-int(c-int(c)<-0.5)

seed=int(sys.argv[1])
e_L=int(sys.argv[2])
n_Dr=int(sys.argv[3])
n_v0=int(sys.argv[4])
e_Dth=int(sys.argv[5])

cla="%d-%d-%d-%d-%d"%(seed,e_L,n_Dr,n_v0,e_Dth)

data_parameter=np.array(pd.read_csv("_bb-20250708-data-%s.txt"%cla,header=None,sep=" ",engine="python",nrows=1))
seed=int(data_parameter[0,0])
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

nB=my_round(L_square/R_inter)
def p(N):
	p1=p0*(1.0-float(N)/Ns)
	if p1>0.0:
		return p1
	else:
		return 0.0
def bc_l(l):
	if l>0.5*L_square:
		return L_square-l
	elif l<-0.5*L_square:
		return L_square+l
	else:
		return l
def bc_b(b):
	if b==-1:
		return nB-1
	elif b==nB:
		return 0
	else:
		return b
def bN(k,y,x):
	if k==0:
		return bc_b(y-1)*nB+bc_b(x-1)
	elif k==1:
		return bc_b(y-1)*nB+bc_b(x*1)
	elif k==2:
		return bc_b(y-1)*nB+bc_b(x+1)
	elif k==3:
		return bc_b(y*1)*nB+bc_b(x-1)
	elif k==4:
		return bc_b(y*1)*nB+bc_b(x*1)
	elif k==5:
		return bc_b(y*1)*nB+bc_b(x+1)
	elif k==6:
		return bc_b(y+1)*nB+bc_b(x-1)
	elif k==7:
		return bc_b(y+1)*nB+bc_b(x*1)
	elif k==8:
		return bc_b(y+1)*nB+bc_b(x+1)
	else:
		return -1

data_movie=np.array(pd.read_csv("_bb-20250708-data-%s.txt"%cla,header=None,sep=" ",engine="python",skiprows=1,skipfooter=2))

fig=plt.figure(figsize=(10,8))
ax=[fig.add_subplot(1,4,(1,3))]
fs=10
lw=0.5
ms=2/L_square
ds=0.02
alpha_c=1.0
alpha_p=1.0
nP=50
lP=[(j+0.5)/nP for j in range(nP)]

ax[0].set_xticks([0,L_square])
ax[0].set_xticklabels(["0.0","%.1f"%L_square],fontsize=fs)
ax[0].set_xlim([0,L_square])
ax[0].set_yticks([0,L_square])
ax[0].set_yticklabels(["0.0","%.1f"%L_square],fontsize=fs)
ax[0].set_ylim([0,L_square])
ax[0].set_aspect(1)

l=0
J_snapshot=0
while l<len(data_movie) and J_snapshot==0:
	tau=int(data_movie[l,0])
	M_now=int(data_movie[l,1])
	g_max=float(data_movie[l,2])
	S1=float(data_movie[l,3])
	l+=1
	if tau==tau_snapshot:
		p_=0.0
		ax[0].set_title("$t=$%.2f, $M=$%d, $g_{max}=$%.2f, $S=$%.2f"%(tau*dt,M_now,g_max,S1),fontsize=fs)
		matrix_p=np.zeros(shape=(nP,nP))
		cB=np.zeros(shape=nB*nB,dtype=int)
		mB=np.zeros(shape=(nB*nB,4*Ns),dtype=int)
		for i in range(M_now):
			x_i=float(data_movie[l+i,1])
			y_i=float(data_movie[l+i,2])
			th_i=float(data_movie[l+i,3])
			N_i=int(data_movie[l+i,4])
			ax[0].plot(x_i,y_i,linestyle="",marker="o",markersize=ms,color="k",alpha=alpha_c)
			ax[0].plot([x_i,x_i+ds*np.cos(th_i)],[y_i,y_i+ds*np.sin(th_i)],linestyle="-",linewidth=lw,color="k",alpha=alpha_c)
			p_+=p(N_i)
			b=int(y_i/R_inter)*nB+int(x_i/R_inter)
			mB[b,cB[b]]=i
			cB[b]+=1
		for j in range(nP):
			for i in range(nP):
				Np=0
				xP=(i+0.5)/nP
				yP=(j+0.5)/nP
				for k in range(9):
					b=bN(k,int(yP/R_inter),int(xP/R_inter))
					for c in range(cB[b]):
						m=mB[b,c]
						x_m=float(data_movie[l+m,1])
						y_m=float(data_movie[l+m,2])
						if np.sqrt(bc_l(xP-x_m)*bc_l(xP-x_m)+bc_l(yP-y_m)*bc_l(yP-y_m))<R_inter:
							Np+=1
				matrix_p[j,i]=p(Np)
		im=ax[0].pcolor(lP,lP,matrix_p,cmap="cool",vmin=0.0,vmax=0.3,alpha=alpha_p,edgecolor="none",shading="nearest")
		cbar=plt.colorbar(im,ax=ax[0],shrink=0.7)
		cbar.ax.tick_params(labelsize=fs)
		cbar.set_label("p(N)",fontsize=fs)
		p_/=M_now
		cbar.ax.plot([0,1],[p_,p_],color="r",linewidth=2)
		J_snapshot=1
	l+=M_now

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
ax[0].plot([],[],label="nP=%d"%nP)
ax[0].legend(loc="upper left",bbox_to_anchor=(1.30,0.00,1,1),ncol=1,handlelength=0,handletextpad=0,fontsize=fs)

fig.savefig("__bb-20250803-t27-%s.pdf"%cla)

