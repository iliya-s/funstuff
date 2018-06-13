import math
import numpy as np
import lanczos as L
import time

#build diagonally dominant Hamiltonian and Direct method function of Hamiltonian
n = 1500
H = np.zeros((n,n))


r = range(n)	#take diagonal elements to be increasing integer values
for i in r:
	H[i,i] = i+1	#take off diagonal elements to be decreasing in order of magnitude
	for j in r[(i+1):]:
		H[i,j] = (10**(i-j+1))
H = (H.T + H)/2
 
def A(v):		#make direct method function
	return np.dot(H,v)

#Lanczos
start_Lanczos = time.time()
 
E = L.lanczos(A,n)

end_Lanczos = time.time()

start_numpy = time.time()

print("Lanczos = ", E,":",end_Lanczos-start_Lanczos, "seconds")

#Numpy
start_numpy = time.time()

E,V = np.linalg.eig(H)
E = np.sort(E)

end_numpy = time.time()

print("Numpy = ", E[0],":",end_numpy-start_numpy,"seconds")
