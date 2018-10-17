# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 08:16:27 2018

@author: xunjie
"""

from __future__ import division
import VfunCPV as Vfun
from scipy.optimize import minimize
import numpy as np
#import time
#
#outpath="v5-result/result"+str(time.time())+".h5"

import sys
istart=0;

if len(sys.argv)==1:
    outpath="v6-result/result-gen"+".h5"
else:
    outpath="v6-result/result-gen"+sys.argv[1]+".h5"
    istart=int(sys.argv[1])
print "outpath=" , outpath


class notBFB():
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

def findmin(c):
    
    x0 = np.random.uniform(-1,1,Vfun.fieldLength)
    # first try to fix phases=0
    fixed_id=[]
    res = minimize(Vfix, x0, args=(c,fixed_id),method='Nelder-Mead', tol=1e-7,
                   options={"maxiter":5000})
    if not(res.success):
        return {"m_L":c,"Vmin":res.fun,"xmin":res.x,"zeroid":[],
                "success":res.success,"flag":0}
    
    #print "success minimized:",res.success
    #print "success minimized:",res.x
    x1=np.copy(res.x);
    fixed_id=np.nonzero(abs(x1)<1e-5)[0]
    
    if len(fixed_id)==0:
        return {"m_L":c,"Vmin":res.fun,"xmin":res.x,"zeroid":[],
                "success":res.success,"flag":1}
    
    res1 = minimize(Vfix, x1, args=(c,fixed_id),method='Nelder-Mead', tol=1e-7,
                   options={"maxiter":1000})
    x2=np.copy(res1.x);
    x2[fixed_id]=0;
    if res1.fun<res.fun:
        return {"m_L":c,"Vmin":res1.fun,"xmin":x2,"zeroid":fixed_id,
                "success":res1.success,"flag":2}
    else:
        return {"m_L":c,"Vmin":res.fun,"xmin":x1,"zeroid":[],
                "success":res.success,"flag":3}

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


    
    
#field_name=["k1","k2","x1","x2","x3","y1","y2","y3",
#            "t2","f1","f2","f3","g1","g2","g3"]

field_name=["k1","k2","x1","x2","x3","y1","y2","y3",
            "","","","","","",""]


#
#
results=[];
n_sample=400000;
seed_list=np.r_[0:n_sample]+istart*n_sample;
for seed in seed_list:
    np.random.seed(seed)
#    print seed;
    c=np.random.uniform(-1,1,Vfun.cLength)
    try:
        checkBFB(c);#if not BFB, it will jump to "excep"
        result=findmin(c)
        result["BFB"]=True       
    except notBFB:
        continue
#   
    result["x8small"]=absmin8(result["xmin"],result["zeroid"]) 
    
    if (result["success"]) and (result["x8small"]>1.0/n_sample):
        print seed, "success";
        result["seed"]=seed
        result["zerofields"] = ''.join(field_name[e] 
                                    for e in result["zeroid"])
        results.append(result)
        
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
f.close()
print "results output:"+outpath
