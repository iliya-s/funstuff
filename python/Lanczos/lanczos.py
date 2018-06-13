import numpy as np
import scipy.linalg

def lanczos(H_func,dimension):
	guess = np.random.rand(dimension)
	q0 = np.zeros(dimension)
	q1 = guess/np.linalg.norm(guess)
	a = []
	b = [] 
	beta = 0.
	for i in range(dimension):
		r = H_func(q1) - beta*q0
		
		alpha = np.dot(q1.T,r)
		r = r - alpha*q1
		beta = np.linalg.norm(r)

		q0 = q1
		q1 = r/beta

		a.append(alpha)
		if i < dimension-1:
			b.append(beta)
	w,v = scipy.linalg.eigh_tridiagonal(a,b)
	w = np.sort(w)
	return w[0]		
