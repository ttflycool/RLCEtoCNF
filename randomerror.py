import numpy as np
from srlce import k,n,w,t
import random

def error(n,w,t):
    e=np.arange(1,n)
    np.random.shuffle(e)
    ee=np.zeros((n+w,),dtype=int)
    for i in range (0,t):
            ee[i]=e[i]
            ee=np.array(ee)
    np.random.shuffle(ee)
    ee=np.array(ee)     
    return ee

rerror=error(n,w,t)