##___bb-20250721_t24.py
## % time python ___bb-20250721_t24.py n_Dr n_v0 e_Dth &

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import glob

def conv(f):
	e=int(f/100)
	for j in range(f%100):
		e*=0.1
	return e;

def arg_th(deg):
	while deg<0:
		deg+=8
	while deg>=8:
		deg-=8
	if deg<=4:
		return deg+5
	elif deg>=5:
		return deg-3
	else:
		return "error"

files=glob.glob("_bb-20250721_t24-data-*.txt")

##switch
n_Dr=int(sys.argv[1])
n_v0=int(sys.argv[2])
e_Dth=int(sys.argv[3])

list_seed=[n for n in range(20)]
Dr=2.0*10.0**(-6.0+n_Dr/4.0)
v0=2.0*10.0**(-4.0+n_v0/4.0)
Dth=conv(e_Dth)
list_factor=[n for n in range(8)]

##array
dth=np.pi*np.linspace(-5,5,11)/4.0
p1=np.zeros(shape=(10,11))
p2=np.zeros(shape=(10,11))
M0=np.zeros(shape=11)
SE=np.zeros(shape=(10,11))
for seed in list_seed:
	for factor in list_factor:
		if "_bb-20250721_t24-data-%d-%d-%d-%d-%d.txt"%(seed,n_Dr,n_v0,e_Dth,factor) in files:
			data=np.array(pd.read_csv("_bb-20250721_t24-data-%d-%d-%d-%d-%d.txt"%(seed,n_Dr,n_v0,e_Dth,factor),header=None,sep=" ",engine="python",skipfooter=1))
			for j in range(len(data)):
				p1[:,arg_th(int(data[j,0]))]+=data[j,2:]
				p2[:,arg_th(int(data[j,0]))]+=data[j,2:]**2
				M0[arg_th(int(data[j,0]))]+=1.0
for j in range(2,10,1):
	p1[:,j]/=M0[j]
	p2[:,j]/=M0[j]
	SE[:,j]=np.sqrt((p2[:,j]-p1[:,j]**2)/M0[j])
p1[:,0]=p1[:,8]
p1[:,1]=p1[:,9]
p1[:,10]=p1[:,2]
M0[0]=M0[8]
M0[1]=M0[9]
M0[10]=M0[2]
SE[:,0]=SE[:,8]
SE[:,1]=SE[:,9]
SE[:,10]=SE[:,2]

##graph
fig=plt.figure(figsize=(20,10))
ax=fig.add_subplot(111)
color_map=plt.get_cmap("gist_rainbow")
fs=30
lw=7
ms=20
tl=12
tw=6
ax.set_title("Dr=%.2e, v0=%.2e, Dth=%.0e"%(Dr,v0,Dth),fontsize=fs)
ax.set_xlabel("$\\Delta\\theta$",fontsize=fs)
ax.set_xticks([-np.pi,-np.pi/2,0,np.pi/2,np.pi])
ax.set_xticklabels(["-$\\pi$","-$\\pi$/2","0","$\\pi$/2","$\\pi$"],fontsize=fs)
ax.set_xlim(-np.pi*9/8,np.pi*9/8)
ax.set_ylabel("$\\langle p\\rangle$",fontsize=fs)
ax.tick_params(axis="both",length=tl,width=tw,labelsize=fs)
for j in range(5):
	ax.vlines(x=dth,ymin=p1[j]-SE[j],ymax=p1[j]+SE[j],color=color_map(j/5),linewidth=lw)
	ax.plot(dth,p1[j],"o-",color=color_map(j/5),markersize=ms,label="%d/q"%(j+1),linewidth=lw)
#ax.legend(loc="best",fontsize=fs)
fig.savefig("__bb-20250721_t24-%d-%d-%d.pdf"%(n_Dr,n_v0,e_Dth))

