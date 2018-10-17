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
    outpath="v8-result/result"+".h5"
else:
    outpath="v8-result/result"+sys.argv[1]+".h5"
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



def findmin(c):
    
    x0 = np.random.uniform(-1,1,Vfun.fieldLength)
    fixed_id=[]
    res = minimize(Vfix, x0, args=(c,fixed_id),method='Nelder-Mead', tol=1e-9,
                   options={"maxiter":5000})
    return {"m_L":c,"Vmin":res.fun,"xmin":res.x,"zeroid":[],
                "success":res.success}
    

def checkBFB(c):
    
    x0 = np.random.uniform(-1,1,Vfun.fieldLength)
    c_nomass=np.copy(c);
    c_nomass[:3]=0;
    fixed_id=[];
    res=minimize(Vfix, x0, args=(c_nomass,fixed_id,-1),method='Nelder-Mead', tol=1e-5,
                   options={"maxiter":1000})
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

cid={'m1q':0,'m2q':1,'m3q':2,
     'L1':3,'L2':4,'L3':5,'L4':6,
     'r1':7,'r2':8,'r3':9,'r4':10,
     'a1':11,'a2':12,'a3':13,'b1':14,'b2':15,'b3':16}

did={'k1':0,'k2':1,
     'x1':2,'x2':3,'x3':4,
     'y1':5,'y2':6,'y3':7}
field_name=["k1","k2","x1","x2","x3","y1","y2","y3",
            "","","","","","",""]

def setsome(c):
#    c[cid['a1']]=0.05*c[cid['a1']];
#    c[cid['a2']]=0.05*c[cid['a2']];
#    c[cid['a3']]=0.05*c[cid['a3']];
#    c[cid['b1']:]=0;#set b1,b2,b3=0
#    c[:]=abs(c[:])
    c[cid['a1']]=0.05*c[cid['a1']];
    c[cid['a2']]=0.05*c[cid['a2']];
    c[cid['a3']]=0.05*c[cid['a3']];

    for positive_label in ['m1q','m2q','m3q','L1','L3','r1','r2','r3']:
        c[cid[positive_label]]=abs(c[cid[positive_label]])   
        
    for label in ['L2','L4','r4','b1','b2','b3']:
        c[cid[label]]=0.0 
    
 
#
lambda_up=4*(np.pi)
lambda_down=-4*(np.pi)

results=[];
n_sample=1000;
seed_list=np.r_[0:n_sample]+istart*n_sample;
for seed in seed_list:
    np.random.seed(seed)
#    print seed;
    c=np.random.uniform(lambda_down,lambda_up,Vfun.cLength)
    setsome(c)
    if c[cid['r3']] <= 2*c[cid['r1']]:
        continue
    
    try:
        checkBFB(c);#if not BFB, it will jump to "excep"
        result=findmin(c)
        result["BFB"]=True       
    except notBFB:
#        print seed,"nBFB"
        continue
#   
   
    if (result["success"]):
        print seed, "success";
        result["seed"]=seed
        result["zerofields"] = ''.join(field_name[e] 
                                    for e in result["zeroid"])
        results.append(result)
        

m_L_data=np.array([results[i]["m_L"] for i in range(len(results))])  
Vmin_data=np.array([results[i]["Vmin"] for i in range(len(results))])  
xmin_data=np.array([results[i]["xmin"] for i in range(len(results))])   
seed_data=np.array([results[i]["seed"] for i in range(len(results))]) 
zerofields_data=  [results[i]["zerofields"] for i in range(len(results))]
#flag_data=  [results[i]["flag"] for i in range(len(results))]
 
     
import h5py
f=h5py.File(outpath, "w")
f.create_dataset("m_L", data=m_L_data)
f.create_dataset("V_min", data=Vmin_data)
f.create_dataset("x_min", data=xmin_data)
f.create_dataset("seed", data=seed_data)
f.create_dataset("zerofields", data=zerofields_data)
#f.create_dataset("flag", data=flag_data)
f.close()
print "results output:"+outpath
