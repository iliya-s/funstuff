#!/usr/bin/env python
import numpy as np
import conjgrad as cj

b = np.array([1,2],dtype = np.float64)
xguess = np.array([0,0],dtype = np.float64)

def Hessian(v):
	A = np.array([[4,1],[1,3]],dtype = np.float64)
	return np.dot(A,v)


x = cj.conjgrad(Hessian,b,xguess)

print(x)
