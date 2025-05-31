## ___bb-20250318-aggregation.py

import numpy as np
import pandas as pd
import sys

seed=int(sys.argv[1])
e_L=int(sys.argv[2])
e_Dr=int(sys.argv[3])
e_v0=int(sys.argv[4])
e_Dth=int(sys.argv[5])

data_parameter=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",nrows=1))
data_movie=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=1,skipfooter=1))
data_diagram=np.array(pd.read_csv("_bb-20250318-data-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),header=None,sep=" ",engine="python",skiprows=len(data_parameter)+len(data_movie)))

np.savetxt("_bb-20250318-data_diagram-%d-%d-%d-%d-%d.txt"%(seed,e_L,e_Dr,e_v0,e_Dth),data_diagram,fmt="%lf")
