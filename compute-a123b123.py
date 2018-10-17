# -*- coding: utf-8 -*-
"""
read inout*.h5
write to result*.h5

@author: xunjie
"""

from __future__ import division
#from IPython import get_ipython
#get_ipython().magic('reset -sf')


import VfunCPV as Vfun
from scipy.optimize import minimize
import numpy as np
import h5py

import sys
istart=0;

if len(sys.argv)==1:
    outpath="scan-a123b123-data/result"+".h5"
    inpath="scan-a123b123-data/input"+".h5"
    min_repeat=3
else:
    outpath="scan-a123b123-data/result"+sys.argv[1]+".h5"
    inpath="scan-a123b123-data/input"+sys.argv[1]+".h5"
    istart=int(sys.argv[1])
    min_repeat=10
print "outpath=" , outpath


class notBFB():
    pass

class unexpected_err():
    pass


def Vfix(x,c,fixed_id,low_bound=-1e5):
    xfix=np.copy(x);
    xfix[fixed_id]=0;
    V=Vfun.V(xfix,c)
    if V>low_bound:
        return V
    else:
        raise notBFB()

phase_id=[8,9,10,11,12,13,14];

  
def findmin_gen(c):
    fixed_id=[]
    x0 = np.random.uniform(-1,1,Vfun.fieldLength)
    res = minimize(Vfix, x0, args=(c,fixed_id),method='Nelder-Mead', tol=1e-7,
                   options={"maxiter":5000})
#    res.x[fixed_id]=0
    return {"m_L":c,"Vmin":res.fun,"xmin":res.x,"zeroid":fixed_id,
                "success":res.success,"flag":0,
                "nit":res.nit, "nfev":res.nfev}
def findmin_success(c):
    while True:
        min=findmin_gen(c);
        if min["success"]:
            return min
            
def findmin_success_n(c,n):
    min0=findmin_success(c)
    for i in range(n):
        min=findmin_success(c)
        if min0["Vmin"]>min["Vmin"]:
            min0=min
    
    return min0;
    

def checkBFB(c):
    
    x0 = np.random.uniform(-1,1,Vfun.fieldLength)
    c_nomass=np.copy(c);
    c_nomass[:3]=0;
    fixed_id=[];
    res=minimize(Vfix, x0, args=(c_nomass,fixed_id,-1),method='Nelder-Mead', tol=1e-7,
                   options={"maxiter":5000})
    #if not BFB, it will raise notBFB() error, 
    #caught by "except" outside
    if res.fun<-1:
        raise notBFB()
    #or if it gets negative V, then BFB is also violated because
    #there is no scales in V ( m^2=0)

    return True          
       
#
def absmin8(xfull,zeroid):
    #find minimal element in 1..8 ,except for zeroid
    xfull2=np.copy(xfull)
    xfull2[zeroid]=1
    x8=xfull2[:8]
    
    return min(abs(x8))

xid={'m1q':0,'m2q':1,'m3q':2,
     'L1':3,'L2':4,'L3':5,'L4':6,
     'r1':7,'r2':8,'r3':9,'r4':10,
     'a1':11,'a2':12,'a3':13,'b1':14,'b2':15,'b3':16}


    
    
#field_name=["k1","k2","x1","x2","x3","y1","y2","y3",
#            "t2","f1","f2","f3","g1","g2","g3"]

field_name=["k1","k2","x1","x2","x3","y1","y2","y3",
            "","","","","","",""]


#
#
fin=h5py.File(inpath, 'r')   
# 'r' means that hdf5 file is open in read-only mode
key=fin.keys()
mLs=np.copy(fin[u'mLs'])
fin.close()

#print np.shape(mLs)

results=[];
m_L_nBFB=[];
np.random.seed(0)

for idx,mL in enumerate(mLs):
#    if idx<1133:
#        continue
    c=np.copy(mL)
    print idx
    try:
        checkBFB(c);#if not BFB, it will jump to "excep"
        result=findmin_success_n(c,min_repeat)
    except notBFB:
        m_L_nBFB.append(c)
        continue
    result["BFB"]=True  
    result["seed"]=idx  
    result["zerofields"] = ''.join(field_name[e] 
                                    for e in result["zeroid"])
    results.append(result)
        
#    for i in range(len(results)):
#        if results[i]['zeroid']==[]:
#            raise unexpected_err()
        
#
#

m_L_data=np.array([results[i]["m_L"] for i in range(len(results))])  
Vmin_data=np.array([results[i]["Vmin"] for i in range(len(results))])  
xmin_data=np.array([results[i]["xmin"] for i in range(len(results))]) 
seed_data=np.array([results[i]["seed"] for i in range(len(results))])   
zerofields_data=  [results[i]["zerofields"] for i in range(len(results))]
flag_data=  [results[i]["flag"] for i in range(len(results))]
 
     
import h5py
f=h5py.File(outpath, "w")
f.create_dataset("m_L", data=m_L_data)
f.create_dataset("V_min", data=Vmin_data)
f.create_dataset("x_min", data=xmin_data)
f.create_dataset("seed", data=seed_data)
f.create_dataset("zerofields", data=zerofields_data)
f.create_dataset("flag", data=flag_data)
f.create_dataset("m_L_nBFB", data=m_L_nBFB)
f.close()
print "results output:"+outpath
