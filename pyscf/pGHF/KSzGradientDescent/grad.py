from pyscf import gto, scf, tools
import numpy as np
import scipy.linalg as lalg
import scipy.optimize as opt
import time
import sys

def contractEri(v2, tdm):
    nso = tdm.shape[0]
    nao = nso // 2

    tdmaa = tdm[:nao, :nao]
    tdmbb = tdm[nao:, nao:]
    tdmab = tdm[:nao, nao:]
    tdmba = tdm[nao:, :nao]

    #the spin indices on the J, K matrices are from the contracted indices
    Jaa = np.einsum('pqrs,sr->pq', v2, tdmaa, dtype = complex, optimize = True)

    Jbb = np.einsum('pqrs,sr->pq', v2, tdmbb, dtype = complex, optimize = True)

    Kaa = np.einsum('psrq,sr->pq', v2, tdmaa, dtype = complex, optimize = True)
    Kbb = np.einsum('psrq,sr->pq', v2, tdmbb, dtype = complex, optimize = True)

    Kab = np.einsum('psrq,sr->pq', v2, tdmba, dtype = complex, optimize = True)
    Kba = np.einsum('psrq,sr->pq', v2, tdmab, dtype = complex, optimize = True)

    G1 = np.zeros((nso, nso), dtype = complex)
    G1[:nao, :nao] = Jaa + Jbb - Kaa
    G1[nao:, nao:] = Jaa + Jbb - Kbb
    G1[:nao, nao:] = - Kba
    G1[nao:, :nao] = - Kab

    return G1

def calcSzSymmetryProjector(nso, sz, nGrid = 8):
    nao = nso // 2
    Wg = []
    Rg = []
    for p in range(nGrid):
        theta = 2 * np.pi * p / nGrid

        w = np.exp(- 1j * theta * sz) / nGrid
        Wg.append(w)

        r = lalg.block_diag(np.exp(1j * theta * 0.5) * np.identity(nao), np.exp(1j * theta * (- 0.5)) * np.identity(nao))
        Rg.append(r)
    return Wg, Rg

def calcEnergy(S, H1, v2, Phi, Wg, Rg):
    #denominator and numerator of energy
    D = 0.0
    N = 0.0

    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.conj().T, S, np.diag(Rg[i]), Phi, dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi, invO, Phi.conj().T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.conj().T, S, np.diag(Rg[i]), Phi.conj(), dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi.conj(), invO, Phi.conj().T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    E = N / D
    return np.real(E)

def calcEnergyAndGradient(S, H1, v2, Phi, Wg, Rg):
    #denominator and numerator of energy
    D = 0.0
    N = 0.0

    #gradient with respect to m,n orbital parameter of the denominator and numerator of energy
    Dmn = np.zeros(Phi.shape, dtype = complex)
    Nmn = np.zeros(Phi.shape, dtype = complex)
    #gradient with respect to m,n* orbital parameter of the denominator and numerator of energy
    Dmn_bar = np.zeros(Phi.shape, dtype = complex)
    Nmn_bar = np.zeros(Phi.shape, dtype = complex)

    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.conj().T, S, np.diag(Rg[i]), Phi, dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to m,n orbital parameter
        detOmn = detO * np.einsum('na,ab,bm,m->mn', invO, Phi.conj().T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOmn_bar = detO * np.einsum('ma,a,ab,bn->mn', S, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi, invO, Phi.conj().T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #drivative of H with respect to m,n orbital parameter
        M1 = H1 + G1

        Amn = np.einsum('na,ab,bm,m->mn', invO, Phi.conj().T, M1, np.diag(Rg[i]), dtype = complex, optimize = True)
        Bmn = np.einsum('na,ab,bc,c,cd,de,ef,fm,m->mn', invO, Phi.conj().T, M1, np.diag(Rg[i]), Phi, invO, Phi.conj().T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        Hmn = Amn - Bmn

        Amn_bar = np.einsum('ma,a,ab,bn->mn', M1, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)
        Bmn_bar = np.einsum('ma,a,ab,bc,cd,de,e,ef,fn->mn', S, np.diag(Rg[i]), Phi, invO, Phi.conj().T, M1, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)
        Hmn_bar = Amn_bar - Bmn_bar

        #averages with symmetry weights
        D += Wg[i] * detO
        Dmn += Wg[i] * detOmn
        Dmn_bar += Wg[i] * detOmn_bar

        N += Wg[i] * detO * H
        Nmn += Wg[i] * (detOmn * H + detO * Hmn)
        Nmn_bar += Wg[i] * (detOmn_bar * H + detO * Hmn_bar)

    #for complex conjugation projection, calculate overlaps with the conjugate of the ket, note that the only nonzero derivatives are with respect to m,n* parameters
    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.conj().T, S, np.diag(Rg[i]), Phi.conj(), dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to m,n orbital parameter
        detOmn_bar = detO * np.einsum('na,ab,bm,m->mn', invO, Phi.conj().T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOmn_bar += detO * np.einsum('ma,a,ab,bn->mn', S, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi.conj(), invO, Phi.conj().T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #drivative of H with respect to m,n orbital parameter
        M1 = H1 + G1

        Amn1_bar = np.einsum('na,ab,bm,m->mn', invO, Phi.conj().T, M1, np.diag(Rg[i]), dtype = complex, optimize = True)
        Amn2_bar = np.einsum('ma,a,ab,bn->mn', M1, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)

        Bmn1_bar = np.einsum('na,ab,bc,c,cd,de,ef,fm,m->mn', invO, Phi.conj().T, M1, np.diag(Rg[i]), Phi.conj(), invO, Phi.conj().T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        Bmn2_bar = np.einsum('ma,a,ab,bc,cd,de,e,ef,fn->mn', S, np.diag(Rg[i]), Phi.conj(), invO, Phi.conj().T, M1, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)

        Hmn_bar = Amn1_bar + Amn2_bar - Bmn1_bar - Bmn2_bar

        #averages with symmetry weights
        D += Wg[i] * detO
        Dmn_bar += Wg[i] * detOmn_bar

        N += Wg[i] * detO * H
        Nmn_bar += Wg[i] * (detOmn_bar * H + detO * Hmn_bar)

    E = N / D
    J = Nmn / D - E * (Dmn / D)
    J_bar = Nmn_bar / D - E * (Dmn_bar / D)

    Jr = np.real(J + J_bar)
    Ji = np.real(1j * (J - J_bar))

    #return E, np.concatenate((np.real(Jr), np.real(Ji)), axis = 0)
    J_vec = np.concatenate((Jr.flatten(), Ji.flatten()), axis = 0)
    return E, J_vec

def calcEnergyAndGradientFiniteDifference(S, H1, v2, Phi, Wg, Rg):
    ds = 1.e-6
    E0 = calcEnergy(S, H1, v2, Phi, Wg, Rg)

    Jr = np.zeros(Phi.shape, dtype = float)
    for a in range(Phi.shape[0]):
        for b in range(Phi.shape[1]):
            Ket = Phi.copy()
            Ket[a, b] += ds
            E = calcEnergy(S, H1, v2, Ket, Wg, Rg)

            Jr[a, b] = (E - E0) / ds

    Ji = np.zeros(Phi.shape, dtype = float)
    for a in range(Phi.shape[0]):
        for b in range(Phi.shape[1]):
            Ket = Phi.copy()
            Ket[a, b] += 1j * ds
            E = calcEnergy(S, H1, v2, Ket, Wg, Rg)

            Ji[a, b] = (E - E0) / ds

    #return E0, np.concatenate((Jr, Ji), axis = 0)
    J_vec = np.concatenate((Jr.flatten(), Ji.flatten()), axis = 0)
    return E, J_vec

#this takes a vector of 2N^2 and turns it into a complex matrix of dimension N
def vec_to_matrix(vec):
    vecr = vec[:len(vec)//2]
    veci = vec[len(vec)//2:]
    N = np.sqrt(len(vecr))
    mat = np.zeros((N, N))
    for m in range(N):
        for n in range(N):
            mat[m, n] = mat[m, n]
    return mat

#this takes a complex matrix of dimension N and turns it into a complex vector of 2N^2
def matrix_to_vec(mat):
    matr = np.real(mat)
    mati = np.imag(mat)
    return np.concatenate((matr.flatten(), mati.flatten()), axis = 0)

#the following four functions are helpers to use scipy's optimization library

def real_to_complex(z): #real vector of length 2n -> complex of length n
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z): #complex vector of length n -> real of length 2n
    return np.concatenate((np.real(z), np.imag(z)))

def fun(params, nelectron, S, H1, v2, Wg, Rg):
    paramVec = real_to_complex(params)
    #wavefunction
    nso = S.shape[0]
    Phi = paramVec.reshape((nso, nelectron))
    return calcEnergy(S, H1, v2, Phi, Wg, Rg)

def jac(params, nelectron, S, H1, v2, Wg, Rg):
    paramVec = real_to_complex(params)
    #wavefunction
    nso = S.shape[0]
    Phi = paramVec.reshape((nso, nelectron))
    J = calcGradient(S, H1, v2, Phi, Wg, Rg)
    return J.flatten()

class DIIS(object):
    counter = 0 #diis does sgd for the first nVec iterations
    nVec = 6

    sizeVector = None
    errorVectors = None

    shapeObject = None
    averagedObjects = None

    def __init__(self, sizeVector, shapeObject, nVec = 6):
        self.sizeVector = sizeVector
        self.errorVectors = []

        self.shapeObject = shapeObject
        self.averagedObjects = []

        self.nVec = nVec
        for i in range(nVec):
            errorVec = np.zeros((sizeVector, 1), dtype = float)
            self.errorVectors.append(errorVec)

            averagedObject = np.zeros(shapeObject, dtype = float)
            self.averagedObjects.append(averagedObject)

    def restart(self):
        self.counter = 0

    def update(self, errVec, avgObj):
        self.counter += 1
        self.errorVectors.pop(0)
        self.errorVectors.append(errVec)
        self.averagedObjects.pop(0)
        self.averagedObjects.append(avgObj)

        returnObj = np.zeros(self.shapeObject, dtype = float)
        if self.counter < self.nVec:
            returnObj = avgObj - errVec
        else:
            b = np.zeros((self.nVec + 1, 1), dtype = float)
            b[self.nVec] = - 1.0

            Bmat = np.zeros((self.nVec + 1, self.nVec + 1), dtype = float)
            for i in range(self.nVec):
                Bmat[self.nVec, i] = - 1.0
                Bmat[i, self.nVec] = - 1.0
            for i in range(self.nVec):
                for j in range(self.nVec):
                    Bmat[i, j] = self.errorVectors[i].conj().T.dot(self.errorVectors[j])

            #w  = lalg.solve(Bmat, b)
            w = lalg.lstsq(Bmat, b)[0]

            w = w[0:self.nVec]
            for i in range(self.nVec):
                returnObj += w[i] * self.averagedObjects[i]

        return returnObj

def pGHF(mol, mo_coeff):
    #basic molecule info and integrals
    sz = float(mol.spin) / 2
    ne = mol.nelectron
    nao = mol.nao
    nso = 2 * mol.nao
    s = mol.intor('int1e_ovlp')
    t = mol.intor('int1e_kin')
    v1 = mol.intor('int1e_nuc')
    v2 = mol.intor('int2e', aosym='s1')

    #spin integrals
    S = lalg.block_diag(s, s)
    h1 = t + v1
    H1 = lalg.block_diag(h1, h1)

    #molecular orbitals
    assert(mo_coeff.shape[0] == nso)
    assert(mo_coeff.shape[1] == nso)
    orbs = mo_coeff + 1j * np.zeros(mo_coeff.shape, dtype = float)

    #Sz symmetry projector
    nGrid = 12
    Wg, Rg = calcSzSymmetryProjector(nso, sz, nGrid)

    #wavefunction
    Phi = orbs[:, 0:ne]
    Phi += 1j * 0.01 * np.random.randn(Phi.shape[0], Phi.shape[1])
    Phi += 0.01 * np.random.randn(Phi.shape[0], Phi.shape[1])

    Eold = 100
    Etol = 1.e-8
    Gtol = 1.e-6
    doPrint = True
    calcStart = time.time()
    diis = DIIS(2 * nso * ne, (2 * nso * ne, ), 6)
    for m in range(50):
        iterStart = time.time()

        #energy
        #E = calcEnergy(S, H1, v2, Phi, Wg, Rg)

        #gradient
        E, G = calcEnergyAndGradient(S, H1, v2, Phi, Wg, Rg)
        Efd, Gfd = calcEnergyAndGradientFiniteDifference(S, H1, v2, Phi, Wg, Rg)
        print("G")
        print(E)
        print(G.T)
        print("G fin-diff")
        print(Efd)
        print(Gfd.T)

        #total energy
        E0 = E + mol.energy_nuc()

        timeEnergyGradient = time.time()

        #calculate update
        params = complex_to_real(Phi.flatten())
        #params = params - 1.0 * G
        params = diis.update(G, params)
        #scipy optimizer
        #params = complex_to_real(Phi.flatten())
        #params = params - dt * Jvec

        #sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'SLSQP', jac = jac, tol = tol)
        #sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'L-BFGS-B', jac = jac, tol = tol)

        #update parameters
        Phi = real_to_complex(params).reshape((nso, ne))
        #if sol.success == True:
        #    Phi = real_to_complex(sol.x).reshape((nso, ne))
        #else:
        #    Phi = real_to_complex(params - dt * Jvec).reshape((nso, ne))

        timeOptimizer = time.time()

        #update error
        Eerror = abs(E - Eold)
        Eold = E
        Gerror = lalg.norm(G)

        #print
        if doPrint == True:
            print(f"-------------------------------- {m} --------------------------------")

            print("Projected values")
            print(f"  Electronic Energy: {E}")
            print(f"  Total Energy: {E0}")
            print(f"  Gradient Norm: {Gerror}")
            print(f"  Energy Error: {Eerror}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

            #print("Scipy Optimizer")
            #print(f"  message: {sol.message}")
            #print(f"  fun: {sol.fun}")
            #print(f"  jac: {lalg.norm(sol.jac)}")
            #print(f"  nit: {sol.nit}")
            #print(f"  Time for Optimizer: {timeOptimizer - timeEnergyGradient}")

            #print("Occupied Orbitals")
            #print(Phi)

        #check for convergence
        if Eerror < Etol or Gerror < Gtol:
            break

    if doPrint == True:
        print(f"\nCalculation Complete")
        print(f"  Total Time: {time.time() - calcStart}")
        print(f"  Total Energy: {E0}")
    return E0, Phi

def readMat(filename, shape, iscomplex):
   if(iscomplex):
     matr = np.zeros(shape)
     mati = np.zeros(shape)
   else:
     mat = np.zeros(shape)
   row = 0
   fileh = open(filename, 'r')
   for line in fileh:
     col = 0
     for coeff in line.split():
       if (iscomplex):
         m = coeff.strip()[1:-1]
         matr[row, col], mati[row, col] = [float(x) for x in m.split(',')]
       else:
         mat[row, col]  = float(coeff)
       col = col + 1
     row = row + 1
   fileh.close()
   if (iscomplex):
     mat = matr + 1j * mati
   return mat

np.set_printoptions(precision=1)
np.set_printoptions(suppress=False)

N = 4
a = 1.4

atomstring = ""
for i in range(N):
    atomstring += f"H 0 0 {i * a};"

mol = gto.Mole()
mol.atom = atomstring
#mol.atom = 'H 0 0 0; H 0 0 1.4; H 0 0 2.8;'
#mol.atom = 'H 0 0 0; H 0 0 1.4'
#mol.atom = 'H 0 0 0'
#mol.atom = 'Li 0 0 0'
mol.basis = 'sto-3g'
#mol.basis = '631g'
#mol.basis = 'ccpvdz'
mol.verbose = 4
mol.spin = 0
mol.unit = 'bohr'
mol.build()

norb = mol.nao
mf = scf.GHF(mol)
dm = mf.get_init_guess()
dm = dm + np.random.rand(2 * norb, 2 * norb) / 100
mf.max_cycle = 100
mf.kernel(dm0 = dm)

S = mf.get_ovlp()
occidx = mf.mo_occ > 0
occOrb = mf.mo_coeff[:, occidx]
#print(occOrb)
fock = mf.get_hcore() + mf.get_veff()
#print(fock)

print("\ninitial solution")
#print(mf.mo_coeff)
E0, mo = pGHF(mol, mf.mo_coeff)

#print("\nkszghf solution")
#mo_coeff = readMat("hf.txt", S.shape, True)
#print(mo_coeff)
#E0, mo = pGHF(mol, mo_coeff)

#print("\n")
#print("Final Result")
#print("Energy")
#print(E0)
#print("Molecular Orbitals")
