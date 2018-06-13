import numpy as np
import Davidson as D
import time

#Build diagonally dominant Hamiltonian and Direct method function of Hamiltonian
n = 1500
H = np.zeros((n,n))


r = range(n)
for i in r:
	H[i,i] = i+1	#take diagonal elements to be increasing integer values
	for j in r[(i+1):]:	#take off diagonal elements to be decreasing in order of magnitude
		H[i,j] = (10**(i-j+1))
H = (H.T + H)/2
 
def A(v):		#make direct method function
	return np.dot(H,v)

#Davidson
start_davidson = time.time()
 
E = D.Davidson(A,np.diag(H),4)

end_davidson = time.time()

start_numpy = time.time()

print("Davidson = ", E,":",end_davidson-start_davidson, "seconds")

#Numpy
start_numpy = time.time()

E,V = np.linalg.eig(H)
E = np.sort(E)

end_numpy = time.time()

print("Numpy = ", E[:4],":",end_numpy - start_numpy, "seconds")
