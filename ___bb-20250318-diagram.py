## ___bb-20250318-diagram.py

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

def ConvertToParameter(f):
    e=int(f/100)
    for j in range(f%100):
        e*=0.1
    return e;

FileInDirectory=glob.glob("_bb-20250318-data_diagram-*.txt")

list_seed=[0,1,2,3,4]
e_L=10
list_Dr=[206,506,105,205,505,104,204]
list_v0=[204,504,103,203,503,102,202]
e_Dth=203

List_Dr=[np.log10(ConvertToParameter(e_Dr)) for e_Dr in list_Dr]
List_v0=[np.log10(ConvertToParameter(e_v0)) for e_v0 in list_v0]
list_M=np.zeros(shape=(len(list_Dr),len(list_v0)))
list_rho=np.zeros(shape=(len(list_Dr),len(list_v0)))
list_g=np.zeros(shape=(len(list_Dr),len(list_v0)))
list_S=np.zeros(shape=(len(list_Dr),len(list_v0)))
list_sampling=np.zeros(shape=(len(list_Dr),len(list_v0)))

for seed in list_seed:
    for e_Dr in list_Dr:
        for e_v0 in list_v0:
            if "_bb-20250318-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth) in FileInDirectory:
                data_diagram=np.array(pd.read_csv("_bb-20250318-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python"))
                list_M[list_Dr.index(e_Dr),list_v0.index(e_v0)]+=float(data_diagram[0,0])
                list_rho[list_Dr.index(e_Dr),list_v0.index(e_v0)]+=float(data_diagram[0,0])/(0.01*e_L**2)
                list_g[list_Dr.index(e_Dr),list_v0.index(e_v0)]+=float(data_diagram[0,-1])-1.0
                list_S[list_Dr.index(e_Dr),list_v0.index(e_v0)]+=float(data_diagram[0,1])
                list_sampling[list_Dr.index(e_Dr),list_v0.index(e_v0)]+=1.0
list_M/=list_sampling
list_rho/=list_sampling
list_g/=list_sampling
list_S/=list_sampling

fig=plt.figure(figsize=(10,10))
ax=[fig.add_subplot(7,7,(1,17)),fig.add_subplot(7,7,(5,21)),fig.add_subplot(7,7,(29,45)),fig.add_subplot(7,7,(33,49))]
fs=8
lw=2
shrink=0.8
shading="nearest" #"flat" "nearest"
for j in range(len(ax)):
    ax[j].set_xlabel("$v_0$",labelpad=0,fontsize=fs)
    ax[j].set_xticks([v0 for v0 in List_v0])
    ax[j].set_xticklabels(["%.0e"%ConvertToParameter(e_v0) for e_v0 in list_v0],fontsize=fs)
    ax[j].set_ylabel("$D_r$",labelpad=0,fontsize=fs)
    ax[j].set_yticks([Dr for Dr in List_Dr])
    ax[j].set_yticklabels(["%.0e"%ConvertToParameter(e_Dr) for e_Dr in list_Dr],fontsize=fs)
    ax[j].set_aspect(1)

#ax[0].set_title(r"$\langle M\rangle$, $L=$%.1f, $D_\theta=$%.0e"%(e_L*0.1,ConvertToParameter(e_Dth)),fontsize=fs)
#image_M=ax[0].pcolormesh(List_v0,List_Dr,list_M,cmap="Reds",vmin=1200.0,vmax=2400.0,shading=shading)
#cbar_M=plt.colorbar(image_M,ax=ax[0],shrink=shrink,location="right",orientation="vertical")
#cbar_M.set_ticks(np.arange(1200,2600,200))
#cbar_M.set_ticklabels(["%.0f"%j for j in np.arange(1200,2600,200)])
#cbar_M.ax.tick_params(labelsize=fs)

ax[0].set_title(r"$\langle\rho\rangle$, $L=$%.1f, $D_\theta=$%.0e"%(e_L*0.1,ConvertToParameter(e_Dth)),fontsize=fs)
image_rho=ax[0].pcolormesh(List_v0,List_Dr,list_rho,cmap="Reds",vmin=1200.0,vmax=2400.0,shading=shading)
cbar_rho=plt.colorbar(image_rho,ax=ax[0],shrink=shrink,location="right",orientation="vertical")
cbar_rho.set_ticks(np.arange(1200,2600,200))
cbar_rho.set_ticklabels(["%.0f"%j for j in np.arange(1200,2600,200)])
cbar_rho.ax.tick_params(labelsize=fs)

ax[1].set_title(r"$g(r_2)-1$, $L=$%.1f, $D_\theta=$%.0e"%(e_L*0.1,ConvertToParameter(e_Dth)),fontsize=fs)
image_g=ax[1].pcolormesh(List_v0,List_Dr,list_g,cmap="Greens",vmin=0.0,vmax=0.6,shading=shading)
cbar_g=plt.colorbar(image_g,ax=ax[1],shrink=shrink,location="right",orientation="vertical")
cbar_g.set_ticks(np.linspace(0,0.6,7))
cbar_g.set_ticklabels(["%.1f"%j for j in np.linspace(0,0.6,7)])
cbar_g.ax.tick_params(labelsize=fs)
if e_L==10:
    ax[1].axhline(y=np.log10(3.78e-05),color="r",linestyle="--",linewidth=lw)

ax[2].set_title(r"$\langle S\rangle$, $L=$%.1f, $D_\theta=$%.0e"%(e_L*0.1,ConvertToParameter(e_Dth)),fontsize=fs)
image_S=ax[2].pcolormesh(List_v0,List_Dr,list_S,cmap="Blues",vmin=0.0,vmax=1.0,shading=shading)
cbar_S=plt.colorbar(image_S,ax=ax[2],shrink=shrink,location="right",orientation="vertical")
cbar_S.set_ticks(np.linspace(0,1,6))
cbar_S.set_ticklabels(["%.1f"%j for j in np.linspace(0,1,6)])
cbar_S.ax.tick_params(labelsize=fs)

ax[3].set_title(r"sampling size, $L=$%.1f, $D_\theta=$%.0e"%(e_L*0.1,ConvertToParameter(e_Dth)),fontsize=fs)
cmap_sampling=plt.get_cmap("tab10")
cmap_sampling.set_under("white")
cmap_sampling.set_over("Black")
image_sampling=ax[3].pcolormesh(List_v0,List_Dr,list_sampling,cmap=cmap_sampling,vmin=0.5,vmax=10.5,shading=shading)
cbar_sampling=plt.colorbar(image_sampling,ax=ax[3],shrink=shrink,location="right",orientation="vertical")
cbar_sampling.set_ticks(np.arange(1,11,1))
cbar_sampling.set_ticklabels(["%.0f"%j for j in np.arange(1,11,1)])
cbar_sampling.ax.tick_params(labelsize=fs)

fig.savefig("__bb-20250318-diagram1-L%d-Dth%d.png"%(e_L,e_Dth))
