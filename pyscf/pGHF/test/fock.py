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
        Phig = Rg[i].dot(Phi)

        #overlap quantities
        O = Phi.conj().T.dot(S).dot(Phig)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #density matrix
        tdm = Phig.dot(invO).dot(Phi.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    E = N / D
    return np.real(E)

def calcFockMatrix(S, H1, v2, Phi, Wg, Rg):
    #deformed density matrix
    dm = Phi.dot(Phi.conj().T)

    #denominator and numerator of energy
    D = 0.0
    N = 0.0
    #gradient with respect to a,b element of deformed denisty matrix of the denominator and numerator of energy
    Dab = np.zeros(S.shape, dtype = complex)
    Nab = np.zeros(S.shape, dtype = complex)

    for i in range(len(Wg)):
        Phig = Rg[i].dot(Phi)

        #overlap quantities
        O = Phi.conj().T.dot(S).dot(Phig)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to a,b dm parameter
        detOab  = detO * np.einsum('ab,bc,cm,nd,d,de,ef,fa->mn', invO, Phi.conj().T, S, S, np.diag(Rg[i]), dm, S, Phi, dtype = complex, optimize = True)
        detOab += detO * np.einsum('ab,bc,cd,de,em,m,nf,fa->mn', invO, Phi.conj().T, S, dm, S, np.diag(Rg[i]), S, Phi, dtype = complex, optimize = True)

        #transition density matrix
        tdm = Phig.dot(invO).dot(Phi.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #drivative of H with respect to a,b dm parameter
        M1 = H1 + G1

        X = np.einsum('ab,bc,cd,de,ef->af', S, Phi, invO, Phi.conj().T, S, dtype = complex, optimize = True)

        A1 = np.einsum('na,ab,bm,m->mn', X, dm, M1, np.diag(Rg[i]))
        A2 = np.einsum('na,a,ab,bm->mn', M1, np.diag(Rg[i]), dm, X)

        B1 = - np.einsum('ab,b,bc,cm,nd,d,de,ef,fa->mn', M1, np.diag(Rg[i]), dm, X, S, np.diag(Rg[i]), dm, X, dm, dtype = complex, optimize = True)
        B2 = - np.einsum('ab,b,bc,cd,de,em,m,nf,fa->mn', M1, np.diag(Rg[i]), dm, X, dm, S, np.diag(Rg[i]), X, dm, dtype = complex, optimize = True)

        Hab = A1 + A2 + B1 + B2

        #averages with symmetry weights
        D += Wg[i] * detO
        Dab += Wg[i] * detOab

        N += Wg[i] * detO * H
        Nab += Wg[i] * (detOab * H + detO * Hab)

    E = N / D
    F = (Nab / D) - E * (Dab / D)

    #separate occupied-virtual spaces
    P = dm.dot(S) #occ projector
    I = np.identity(P.shape[0], dtype = complex)
    Q = I - P #virt projector
    Fstandard = H1 + contractEri(v2, dm)
    F = F + P.conj().T.dot(Fstandard).dot(P) + Q.conj().T.dot(Fstandard).dot(Q)
    return F

def calcEnergyFromDensity(S, H1, v2, Phi, dm, Wg, Rg):
    #denominator and numerator of energy
    D = 0.0
    N = 0.0

    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('xa,ab,bc,cd,d,de,ef,fy->xy', Phi.conj().T, S, dm, S, np.diag(Rg[i]), dm, S, Phi, dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #density matrix
        tdm = np.einsum('x,xa,ab,bc,cd,de,ef,fy->xy', np.diag(Rg[i]), dm, S, Phi, invO, Phi.conj().T, S, dm, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    E = N / D
    return np.real(E)

def calcFockMatrixFiniteDifference(S, H1, v2, Phi, Wg, Rg):
    ds = 1.e-6

    dm = Phi.dot(Phi.conj().T)
    E0 = calcEnergyFromDensity(S, H1, v2, Phi, dm, Wg, Rg)

    Fr = np.zeros(dm.shape, dtype = complex)
    for a in range(dm.shape[0]):
        for b in range(dm.shape[1]):
            dm1 = dm.copy()
            dm1[a, b] += ds
            E = calcEnergyFromDensity(S, H1, v2, Phi, dm1, Wg, Rg)

            Fr[a, b] = (E - E0) / ds

    Fi = np.zeros(dm.shape, dtype = complex)
    for a in range(dm.shape[0]):
        for b in range(dm.shape[1]):
            dm1 = dm.copy()
            dm1[a, b] += 1j * ds
            E = calcEnergyFromDensity(S, H1, v2, Phi, dm1, Wg, Rg)

            Fi[a, b] = (E - E0) / ds

    return Fr + 1j * Fi

def geig(A, B):
    w, v = lalg.eigh(B)
    idx = []
    for i in range(w.size):
        if abs(w[i]) > 1.e-8:
            idx.append(i)
    X = np.zeros((B.shape[0], len(idx)), dtype = complex)
    for i in range(len(idx)):
        X[:, i] = v[:, idx[i]] / np.sqrt(w[idx[i]])

    Ap = X.conj().T.dot(A).dot(X)
    Bp = X.conj().T.dot(B).dot(X)
    w, v = lalg.eig(Ap)
    v = X.dot(v)
    return w, v

def geigh(A, B):
    w, v = lalg.eigh(B)
    idx = []
    for i in range(w.size):
        if abs(w[i]) > 1.e-8:
            idx.append(i)
    X = np.zeros((B.shape[0], len(idx)), dtype = complex)
    for i in range(len(idx)):
        X[:, i] = v[:, idx[i]] / np.sqrt(w[idx[i]])

    Ap = X.conj().T.dot(A).dot(X)
    Bp = X.conj().T.dot(B).dot(X)
    w, v = lalg.eigh(Ap)
    v = X.dot(v)
    return w, v

def randOrthoMat(dim, S = None):
    randMat = np.random.randn(dim, dim) + 1j * np.random.randn(dim, dim)
    randMat = randMat + randMat.conj().T
    w, v = lalg.eigh(randMat, S)
    return v

def modifiedGramSchmidt(V, S = None):
    if S is None:
        S = np.identity(V.shape[0], dtype = V.dtype)
    U = np.zeros(V.shape, dtype = V.dtype)
    U[:, 0] = V[:, 0] / np.sqrt( V[:, 0].conj().T.dot(S).dot(V[:, 0]) )
    for i in range(1, V.shape[1]):
        U[:, i] = V[:, i]
        for j in range(i):
            U[:, i] = U[:, i] - ( U[:, j].conj().T.dot(S).dot(U[:, i]) ) / ( U[:, j].conj().T.dot(S).dot(U[:, j]) ) * U[:, j]
        U[:, i] = U[:, i] / np.sqrt( U[:, i].conj().T.dot(S).dot(U[:, i]) )
    return U

class DIIS(object):
    counter = 0 #diis does nothing for the first nVec iterations
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
            errorVec = np.zeros((sizeVector, 1), dtype = complex)
            self.errorVectors.append(errorVec)

            averagedObject = np.zeros(shapeObject, dtype = complex)
            self.averagedObjects.append(averagedObject)

    def restart(self):
        self.counter = 0

    def update(self, errVec, avgObj):
        self.counter += 1
        self.errorVectors.pop(0)
        self.errorVectors.append(errVec)
        self.averagedObjects.pop(0)
        self.averagedObjects.append(avgObj)

        returnObj = np.zeros(self.shapeObject, dtype = complex)
        if self.counter < self.nVec:
            returnObj = avgObj
        else:
            b = np.zeros((self.nVec + 1, 1), dtype = complex)
            b[self.nVec] = - 1.0

            Bmat = np.zeros((self.nVec + 1, self.nVec + 1), dtype = complex)
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

class AMSgrad(object):
    stepsize = 0.001
    decay1 = 0.1
    decay2 = 0.001
    mom1 = None
    mom2 = None

    def __init__(self, stepsize, decay1, decay2, nParams):
        self.stepsize = stepsize
        self.decay1 = decay1
        self.decay2 = decay2
        self.mom1 = np.zeros((nParams, 1), dtype = float)
        self.mom2 = np.zeros((nParams, 1), dtype = float)

    def update(self, gradient):
        assert(gradient.shape[0] == self.mom1.shape[0])

        delta = np.zeros(gradient.shape)
        for i in range(gradient.shape[0]):
            self.mom1[i] = self.decay1 * gradient[i] + (1.0 - self.decay1) * self.mom1[i]
            self.mom2[i] = max(self.mom2[i], self.decay2 * gradient[i] * gradient[i] + (1.0 - self.decay2) * self.mom2[i])

            delta[i] = - self.stepsize * self.mom1[i] / (np.sqrt(self.mom2[i]) + 1.e-8)

        return delta

def pGHF(mol, mo_coeff, ovlp_mat = None):
    #basic molecule info and integrals
    sz = float(mol.spin) / 2
    ne = mol.nelectron
    nao = mol.nao
    nso = 2 * mol.nao
    s = mol.intor('int1e_ovlp')
    t = mol.intor('int1e_kin')
    v1 = mol.intor('int1e_nuc')
    v2 = mol.intor('int2e', aosym='s1')
    assert(mo_coeff.shape[0] == nso)
    assert(mo_coeff.shape[1] == nso)

    #spin integrals
    S = lalg.block_diag(s, s)
    h1 = t + v1
    H1 = lalg.block_diag(h1, h1)

    #molecular orbitals
    orbs = mo_coeff + 1j * np.zeros(mo_coeff.shape, dtype = float)

    #mixes homo and lumo molecular orbitals to break complex conjugation symmetry
    randOccVirMat = randOrthoMat(2)
    randU = np.identity(nso, dtype = complex)
    randU[(ne - 1):(ne + 1), (ne - 1):(ne + 1)] = randOccVirMat
    #orbs = orbs.dot(randU)

    #Sz symmetry projector
    nGrid = 10
    Wg, Rg = calcSzSymmetryProjector(nso, sz, nGrid)

    Eold = 100
    Etol = 1.e-8
    Gtol = 1.e-6
    doPrint = True
    calcStart = time.time()
    diis = DIIS(2 * (nso - ne) * ne, (nso, nso), 6)
    for m in range(50):
        iterStart = time.time()

        #current wavefunction
        Phi = orbs[:, 0:ne]

        #energy
        E = calcEnergy(S, H1, v2, Phi, Wg, Rg)

        #total energy
        E0 = E + mol.energy_nuc()

        #fock matrix
        F = calcFockMatrix(S, H1, v2, Phi, Wg, Rg)
        Fmo = orbs.conj().T.dot(F).dot(orbs)

        #F_fd = calcFockMatrixFiniteDifference(S, H1, v2, Phi, Wg, Rg)
        #Fmo_fd = orbs.conj().T.dot(F_fd).dot(orbs)

        #gradient
        grad = Fmo[0:ne, ne:nso].flatten()
        gradReal = np.concatenate((np.real(grad), np.imag(grad)), axis = 0)

        timeEnergyGradient = time.time()

        #diis
        F = diis.update(gradReal, F)

        #diagonalize fock matrix
        w, v = geigh(F, S)

        #update orbitals if new orbitals lower the energy
        Phip = v[:, 0:ne]
        Ep = calcEnergy(S, H1, v2, Phip, Wg, Rg)
        if Ep < E:
            orbs = v

            #update error
            Eerror = abs(E - Eold)
            Eold = E

            Gerror = lalg.norm(gradReal)

        #print
        if doPrint == True:
            print(f"-------------------------------- {m} --------------------------------")

            print(f"  Electronic Energy: {E}")
            print(f"  Total Energy: {E0}")
            print(f"  Gradient Norm: {lalg.norm(gradReal)}")
            print(f"  Energy Error: {Eerror}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

        #check for convergence
        if Eerror < Etol or Gerror < Gtol:
            break

    if doPrint == True:
        print(f"\nCalculation Complete")
        print(f"  Total Time: {time.time() - calcStart}")
        print(f"  Total Energy: {E0}")
    return E0, orbs

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

np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

N = 3
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
mol.spin = 1
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
print(fock)
print("\n")

E0, mo = pGHF(mol, mf.mo_coeff)

#mo_coeff = readMat("hf.txt", S.shape, True)
#print(mo_coeff)
#E0, mo = pGHF(mol, mo_coeff)

#print("\n")
#print("Final Result")
#print("Energy")
#print(E0)
#print("Molecular Orbitals")
#print(mo)
