#!/usr/bin/env python

import math
import numpy as np
import time

n = 1200 			#dimension of matrix
tol = 1e-7 		#convergence criterion
sparsity = 0.001 	#how diagonally dominant your matrix is

#Build diagonally dominant matrix
H = np.zeros((n,n))

for i in range(n): #Take diagonal values to be increasing integers
	H[i,i] = i + 1

H = H + sparsity*np.random.rand(n,n)
H = (H.T + H)/2


#Initialize
a = 12 			#number of vectors in initial sample space
neig = 4		#number of eigenvalues to solve for
t = np.eye(n,a) 	#set of test unit vectors in initial sample space
I = np.identity(n) 	#identity matrix with dimension of matrix in question
V = np.zeros((n,n)) 	#array to store sample space

for i in range(a):	#input test vectors into sample space matrix
	V[:,i] = t[:,i]


#Davidson

start_davidson = time.time()

theta_old = np.zeros(neig)	#initialize old and new eigenvalue guesses, "Theta"
theta_new = np.ones(neig)
count = 1			#keep track of number of iterations
while np.linalg.norm(theta_old-theta_new) > tol:
	theta_old = theta_new	#step theta
	V,R = np.linalg.qr(V)	#use python's QR decomp. to ensure sample space orthogonality
	HV = np.dot(H,V[:,:(a*count)])
	VHV = np.dot(V[:,:(a*count)].T,HV) #build matrix in subspace
	theta,s = np.linalg.eig(VHV)	#diagonalize
	index = np.argsort(theta)	#sort eigenvalues and eigenvectors
	theta = theta[index]
	s = s[:,index]
	for i in range(a):		#loop through test vectors
		test = np.dot(V[:,:(a*count)],s[:,i])	#change basis of eigenvectors into basis of original matrix
		b = theta[i]*I
		r = np.dot(H,test) - np.dot(b,test)	#calculate residue vector
		q = -(1/H[i,i] - 1/theta[i])*r		#calculate correction vector
		V[:,(i+(a*count))] = q			#add correction vectors to subspace
	theta_new = theta[:neig]			#update guesses to eigenvalues
	count = count + 1

end_davidson = time.time()


print "Davidson = ", theta_new[:neig],";",end_davidson-start_davidson, "seconds"

#Numpy

start_numpy = time.time()

E,V = np.linalg.eig(H)
E = np.sort(E)

end_numpy = time.time()

print "Numpy = ", E[:neig],";",end_numpy - start_numpy, "seconds"
