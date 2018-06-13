import numpy as np

def jacobi(H_func,diag_H,b):
	d = len(diag_H)
	xguess = np.random.rand(d)
	while True:
		R_x = H_func(xguess) - diag_H*xguess
		A = b - R_x
		x = A/diag_H
		if np.linalg.norm(x-xguess) < 1e-8:
			break
		xguess = x
	return x
