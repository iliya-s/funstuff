#!/usr/bin/env python
import numpy as np
import jacobi as j

b = np.array([11,13],dtype = np.float64)
xguess = np.array([1,1],dtype = np.float64)
A = np.array([[2,1],[5,7]],dtype = np.float64)

def Hessian(v):
	return np.dot(A,v)


x = j.jacobi(Hessian,np.diag(A),b)

print(x)
