## ___bb-20250720-L_dependency.py
## % time python ___bb-20250720-L_dependency.py n_Dr e_Dth &

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
import sys

def conv(f):
	e=int(f/100)
	for j in range(f%100):
		e*=0.1
	return e

def XIn1(x):
	while x<0.0:
		x+=1.0
	while x>=1.0:
		x-=1.0
	return x

files=glob.glob("_bb-20250711-data_diagram-*.txt")

#range_v0=[n for n in range(2,5,1)]
range_v0=[2,3,3.5,4]

list_seed=[n for n in range(10)]
list_L=[10,15,20]
n_Dr=int(sys.argv[1])
list_v0=[np.log10(2.0)+(-4.0+n/4.0) for n in range_v0]
e_Dth=int(sys.argv[2])

list_S1=np.zeros(shape=(len(list_L),len(list_v0)))
list_S2=np.zeros(shape=(len(list_L),len(list_v0)))
list_g1=np.zeros(shape=(len(list_L),len(list_v0)))
list_g2=np.zeros(shape=(len(list_L),len(list_v0)))
list_sampling=np.zeros(shape=(len(list_L),len(list_v0)))
for seed in list_seed:
	for e_L in list_L:
		for n_v0 in range_v0:
			J35=0
			if (n_v0!=3.5)and("_bb-20250711-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,n_Dr,n_v0,e_Dth) in files):
				data=np.array(pd.read_csv("_bb-20250711-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,n_Dr,n_v0,e_Dth),header=None,sep=" ",engine="python"))
				J35=1
			if (n_v0==3.5)and("_bb-20250711-data_diagram-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,n_Dr,n_v0,e_Dth) in files):
				data=np.array(pd.read_csv("_bb-20250711-data_diagram-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,n_Dr,n_v0,e_Dth),header=None,sep=" ",engine="python"))
				J35=1
			if J35==1:
				list_S1[list_L.index(e_L),range_v0.index(n_v0)]+=float(data[0,1])
				list_S2[list_L.index(e_L),range_v0.index(n_v0)]+=float(data[0,1])**2
				list_g1[list_L.index(e_L),range_v0.index(n_v0)]+=float(data[0,-1])-1.0
				list_g2[list_L.index(e_L),range_v0.index(n_v0)]+=(float(data[0,-1])-1.0)**2
				list_sampling[list_L.index(e_L),range_v0.index(n_v0)]+=1.0
list_S1/=list_sampling
list_S2/=list_sampling
list_S_SE=np.sqrt((list_S2-list_S1**2)/list_sampling)
list_g1/=list_sampling
list_g2/=list_sampling
list_g_SE=np.sqrt((list_g2-list_g1**2)/list_sampling)

fs=10
lw=2
fig=plt.figure(figsize=(10,6))
ax=[fig.add_subplot(8,5,(1,7)),fig.add_subplot(8,5,(4,10)),fig.add_subplot(8,5,(16,22)),fig.add_subplot(8,5,(19,25)),fig.add_subplot(8,5,(31,37))]
cmap=plt.get_cmap("tab10")
alpha=1.0

for j in range(len(ax)):
	ax[j].set_xlabel("$v_0$",labelpad=0,fontsize=fs)
	ax[j].set_xticks(list_v0)
	ax[j].set_xticklabels(["%.3e"%(10.0**v0) for v0 in list_v0],rotation=90,fontsize=fs)
	if j<2:
		ax[j].set_title("$D_r=%.3e,D_\\theta=%.0e$"%(2.0*10.0**(-6.0+n_Dr/4.0),conv(e_Dth)),fontsize=fs)
	else:
		ax[j].set_ylabel("$L$",labelpad=0,fontsize=fs)
		ax[j].set_yticks(np.array(list_L)*0.1)
		ax[j].set_yticklabels(["%.1f"%(e_L*0.1) for e_L in list_L],fontsize=fs)
		ax[j].set_ylim([0,4])
		ax[j].plot(list_v0,np.zeros(shape=len(list_v0))-1,"k")

ax[0].set_ylabel("$g(r_2)-1$",labelpad=0,fontsize=fs)
ax[0].tick_params(axis="y",labelsize=fs)
ax[1].set_ylabel("$\\langle S\\rangle$",labelpad=0,fontsize=fs)
ax[1].tick_params(axis="y",labelsize=fs)

for n_L in range(len(list_L)):
	ax[0].plot(list_v0,list_g1[n_L],color=cmap(XIn1(n_L*0.1+0.05)),linestyle="-",linewidth=lw,label="$L=%.1f$"%(list_L[n_L]*0.1))
	ax[0].vlines(x=list_v0,ymin=list_g1[n_L]-list_g_SE[n_L],ymax=list_g1[n_L]+list_g_SE[n_L],color=cmap(XIn1(n_L*0.1+0.05)),linestyle="-",linewidth=lw)
	ax[1].plot(list_v0,list_S1[n_L],color=cmap(XIn1(n_L*0.1+0.05)),linestyle="-",linewidth=lw,label="$L=%.1f$"%(list_L[n_L]*0.1))
	ax[1].vlines(x=list_v0,ymin=list_S1[n_L]-list_S_SE[n_L],ymax=list_S1[n_L]+list_S_SE[n_L],color=cmap(XIn1(n_L*0.1+0.05)),linestyle="-",linewidth=lw)
	[ax[2].text(list_v0[range_v0.index(n_v0)],0.1*list_L[n_L],"%.1e"%list_g_SE[n_L,range_v0.index(n_v0)],fontsize=6) for n_v0 in range_v0]
	[ax[3].text(list_v0[range_v0.index(n_v0)],0.1*list_L[n_L],"%.1e"%list_S_SE[n_L,range_v0.index(n_v0)],fontsize=6) for n_v0 in range_v0]
	[ax[4].text(list_v0[range_v0.index(n_v0)],0.1*list_L[n_L],"%.0f"%list_sampling[n_L,range_v0.index(n_v0)],fontsize=fs) for n_v0 in range_v0]

ax[0].legend(loc="upper right",fontsize=fs)
ax[1].legend(loc="upper left",fontsize=fs)

#print("g:n=3  :%.3f\ng:n=3.5:%.3f\nS:n=3  :%.3f\nS:n=3.5:%.3f"%(list_g1[0,range_v0.index(3)],list_g1[0,range_v0.index(3.5)],list_S1[0,range_v0.index(3)],list_S1[0,range_v0.index(3.5)]))

fig.savefig("__bb-20250720-L_dependency-%d-%d.pdf"%(n_Dr,e_Dth))

