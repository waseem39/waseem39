# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 22:38:05 2024

@author: Wasim
"""

from sympy import *


# Jocabian

def Joc(F,Var,N):
    A=zeros(N,N)
    for i in range(0,N):
        for j in range(0,N):
            A[i,j]=diff(F[i],Var[j])
    
    return A


# Multilinear form B

def MulB(F,Var,Var1,N):
    B=zeros(N,1)

    for i in range(0,N):
        s=0
        for j in range(0,N):
            for k in range(0,N):
                dd=diff(F[i],Var[j],Var[k])
                for t in range(0,N):
                    dd=dd.subs({Var[t]:0})
                
                
                s=s+(dd*Var[j]*Var1[k])
        
        B[i]=s
        
    return B


# Multilinear form C

def MulC(F,Var,Var1,Var2,N):
    C=zeros(N,1)

    for i in range(0,N):
        s=0
        for j in range(0,N):
            for k in range(0,N):
                for l in range(0,N):    
                    dd=diff(F[i],Var[j],Var[k],Var[l])
                    for t in range(0,N):
                        dd=dd.subs({Var[t]:0})
                        
                    
                    s=s+(dd*Var[j]*Var1[k]*Var2[l])
        
        C[i]=s
        
    return C

