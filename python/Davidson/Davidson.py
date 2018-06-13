import math
import numpy as np

def Davidson(H_func,H_diag,neig):	#input Direct Method Hamiltonian, Diagonal of Hamiltonian, and number of eigenvalues desired	
	a = 8				#number of vectors in initial sample space
	tol = 1e-8
	n = len(H_diag)
	t = np.eye(n,a)			#set of test unit vectors in initial sample space
	V = np.zeros((n,a))		#array to store sample space

	for i in range(a):		#input test vectors into sample space matrix
		V[:,i] = t[:,i]


	theta_old = np.zeros(neig)	#initialize old and new eigenvalue guesses, "Theta"
	theta_new = np.ones(neig)

	count = 1			#keep track of number of iterations
	while np.linalg.norm(theta_old-theta_new) > tol:
		theta_old = theta_new	#step theta
		V,R = np.linalg.qr(V)	#use python's QR decomp. to ensure sample space orthogonality
		HV = np.zeros((n,a*count))
		for i in range(a*count):
			HV[:,i] = H_func(V[:,i])
		VHV = np.dot(V[:,:(a*count)].T,HV) #build matrix in subspace
		theta,s = np.linalg.eig(VHV)	#diagonalize
		index = np.argsort(theta)	#sort eigenvalues and eigenvectors
		theta = theta[index]
		s = s[:,index]
		V = np.c_[V,np.zeros((n,a))]	#grow sample space matrix
		for i in range(a):		#loop through test vectors
			test = np.dot(V[:,:(a*count)],s[:,i])	#change basis of eigenvectors into basis of original matrix
			r = H_func(test) - theta[i]*test	#calculate residue vector
			q = -(1/H_diag[i] - 1/theta[i])*r	#calculate correction vector
			V[:,(i+(a*count))] = q			#add correction vectors to subspace
		theta_new = theta[:neig]			#update guesses to eigenvalues
		count = count + 1
	return theta_new
