## ___bb-20250711-aggregation.py
## % time python ___bb-20250711-aggregation.py seed e_L n_Dr n_v0 e_Dth &
## % for c in {0..8};do for d in {0..8};do python ___bb-20250711-aggregation.py 0 10 $c $d 203;done;done &
## $ for a in {0..4};do for b in 10;do for c in {0..8};do for d in {0..8};do for e in 203;do bsub -q "jobLimit52" "python ___bb-20250711-aggregation.py $a $b $c $d $e";done;done;done;done;done

import numpy as np
import pandas as pd
import sys

seed=int(sys.argv[1])
e_L=int(sys.argv[2])
e_Dr=int(sys.argv[3])
e_v0=float(sys.argv[4])
e_Dth=int(sys.argv[5])

if e_v0!=3.5:
	data_parameter=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
	data_movie=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=2))
	data_diagram=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=len(data_parameter)+len(data_movie),skipfooter=1))
	np.savetxt("_bb-20250711-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),data_diagram,fmt="%lf")
elif e_v0==3.5:
	data_parameter=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
	data_movie=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=2))
	data_diagram=np.array(pd.read_csv("_bb-20250708-data-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=len(data_parameter)+len(data_movie),skipfooter=1))
	np.savetxt("_bb-20250711-data_diagram-%d-%d-%d-%.1f-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),data_diagram,fmt="%lf")


