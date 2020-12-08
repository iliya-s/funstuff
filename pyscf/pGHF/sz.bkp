from pyscf import gto, scf, tools
import numpy as np
import scipy.linalg as lalg
import scipy.optimize as opt
import time
import sys

def calcOverlap(S, Bra, Ket):
    ovlpMat = Bra.conj().T.dot(S).dot(Ket)
    return lalg.det(ovlpMat)

def calcDensityMatrix(S, Bra, Ket):
    ovlpMat = Bra.conj().T.dot(S).dot(Ket)
    ovlpInv = lalg.inv(ovlpMat)
    return Ket.dot(ovlpInv).dot(Bra.conj().T)

def calcHamiltonianMatrixElement(H1, v2, dm1):
    nso = H1.shape[0]
    nao = nso // 2

    dm1aa = dm1[:nao, :nao]
    dm1bb = dm1[nao:, nao:]
    dm1ab = dm1[:nao, nao:]
    dm1ba = dm1[nao:, :nao]

    #the spin indices on the J, K matrices are from the contracted indices
    Jaa = np.einsum('pqrs,sr->pq', v2, dm1aa, dtype = complex)

    Jbb = np.einsum('pqrs,sr->pq', v2, dm1bb, dtype = complex)

    Kaa = np.einsum('psrq,sr->pq', v2, dm1aa, dtype = complex)
    Kbb = np.einsum('psrq,sr->pq', v2, dm1bb, dtype = complex)

    Kab = np.einsum('psrq,sr->pq', v2, dm1ba, dtype = complex)
    Kba = np.einsum('psrq,sr->pq', v2, dm1ab, dtype = complex)

    G1 = np.zeros((nso, nso), dtype = complex)
    G1[:nao, :nao] = Jaa + Jbb - Kaa
    G1[nao:, nao:] = Jaa + Jbb - Kbb
    G1[:nao, nao:] = - Kba
    G1[nao:, :nao] = - Kab

    H = H1 + 0.5 * G1
    return np.einsum('pq,qp->', H, dm1, dtype = complex)

def calcSzSymmetryProjector(nao, sz, nGrid = 8):
    Wg = []
    Rg = []
    for p in range(nGrid):
        theta = 2 * np.pi * p / nGrid

        w = np.exp(- 1j * theta * sz) / nGrid
        Wg.append(w)

        r = lalg.block_diag(np.exp(1j * theta * 0.5) * np.identity(nao), np.exp(1j * theta * (- 0.5)) * np.identity(nao))
        Rg.append(r)
    return Wg, Rg

def calcEnergyNumeratorDenominator(S, H1, v2, Bra, Ket, Wg, Rg):
    #integrate
    D = 0.0
    N = 0.0
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket)

        #calculate quantities
        Og = calcOverlap(S, Bra, Ketg)

        dmg = calcDensityMatrix(S, Bra, Ketg)
        Eg = calcHamiltonianMatrixElement(H1, v2, dmg)
        Eg = Og * Eg

        #averages
        D += Wg[i] * Og
        N += Wg[i] * Eg
    return N, D

def generateSingleExcitations(ne, nso):
    hf = [x for x in range(ne)]
    excitations = []
    indices = []
    for i in range(ne):
        for a in range(ne, nso):
            #excitation
            excitations.append([i, a])

            #indices
            orbs = hf.copy()
            orbs[i] = a
            indices.append(orbs)
    return excitations, indices

#note this returns an nso by nso matrix but only the top right block is non-zero
def calcGradientEnergyNumeratorDenominator(S, H1, v2, nelectron, mo_coeff, Wg, Rg):
    Npq = np.zeros(S.shape, dtype = complex)
    Dpq = np.zeros(S.shape, dtype = complex)

    #hf det
    Psi = mo_coeff[:, 0:nelectron]

    #single excitations
    excitations, indices = generateSingleExcitations(nelectron, S.shape[0])
    for n in range(len(excitations)):
        Psi_n = mo_coeff[:, indices[n]]
        N, D = calcEnergyNumeratorDenominator(S, H1, v2, Psi_n, Psi, Wg, Rg)

        Npq[tuple(excitations[n])] = N
        Dpq[tuple(excitations[n])] = D
    return Npq, Dpq

def modifiedGramSchmidt(V, S = None):
    U = np.zeros(V.shape)
    if S is None:
        S = np.identity(U.shape[0])

    for i in range(V.shape[1]):
        v = V[:, i]
        for j in range(i):
            u = U[:, j]
            ovlp = v.T.dot(S).dot(u)
            norm2 = u.T.dot(S).dot(u)
            v = v - u * ovlp / norm2

        norm = v.T.dot(S).dot(v)
        U[:, i] = v / np.sqrt(norm)
    return U


#the following four functions are helpers to use scipy's optimization library

def real_to_complex(z): #real vector of length 2n -> complex of length n
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z): #complex vector of length n -> real of length 2n
    return np.concatenate((np.real(z), np.imag(z)))

def fun(params, S, H1, v2, nelectron, mo_coeff, Wg, Rg):
    nso = S.shape[0]
    nao = nso // 2

    #build kappa matrix from params
    vec = real_to_complex(params)
    paramMat = vec.reshape((nelectron, nso - nelectron))
    K = np.zeros(((nso, nso)), dtype = complex)
    K[0:nelectron, nelectron:nso] = paramMat
    K = K - K.conj().T

    #orbital rotation
    U = lalg.expm(- K)
    orbs = mo_coeff.dot(U)

    #wavefunction
    Psi = orbs[:, 0:nelectron]

    #energy terms
    N, D = calcEnergyNumeratorDenominator(S, H1, v2, Psi, Psi, Wg, Rg)
    N = np.real(N)
    D = np.real(D)
    return N / D

def jac(params, S, H1, v2, nelectron, mo_coeff, Wg, Rg):
    nso = S.shape[0]
    nao = nso // 2

    #build kappa matrix from params
    vec = real_to_complex(params)
    paramMat = vec.reshape((nelectron, nso - nelectron))
    K = np.zeros(((nso, nso)), dtype = complex)
    K[0:nelectron, nelectron:nso] = paramMat
    K = K - K.conj().T

    #orbital rotation
    U = lalg.expm(- K)
    orbs = mo_coeff.dot(U)

    #wavefunction
    Psi = orbs[:, 0:nelectron]

    #energy terms
    N, D = calcEnergyNumeratorDenominator(S, H1, v2, Psi, Psi, Wg, Rg)
    N = np.real(N)
    D = np.real(D)

    #gradient terms
    Npq, Dpq = calcGradientEnergyNumeratorDenominator(S, H1, v2, nelectron, orbs, Wg, Rg)

    G = (Npq / D) - (N / D) * (Dpq / D)
    G = G[0:nelectron, nelectron:nso]
    G = G.flatten()
    return complex_to_real(G)

def pGHF(mol, mo_coeff = None):
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

    if mo_coeff is None:
        orbs = np.random.randn(nso, nso)
        orbs = modifiedGramSchmidt(orbs, S)
    else:
        assert(mo_coeff.shape[0] == nso)
        assert(mo_coeff.shape[1] == nso)
        orbs = mo_coeff

    #slight amount of noise helps optimization
    orbs += np.random.randn(nso, nso) / 100

    #Sz symmetry projector
    nGrid = 10
    Wg, Rg = calcSzSymmetryProjector(nao, sz, nGrid)

    Eold = 100
    dt = 0.1
    tol = 1.e-8
    doPrint = True
    calcStart = time.time()
    for m in range(50):
        iterStart = time.time()

        #current wavefunction
        Psi = orbs[:, 0:ne]

        #energy terms
        N, D = calcEnergyNumeratorDenominator(S, H1, v2, Psi, Psi, Wg, Rg)
        N = np.real(N)
        D = np.real(D)

        #gradient terms
        Npq, Dpq = calcGradientEnergyNumeratorDenominator(S, H1, v2, ne, orbs, Wg, Rg)

        #electronic and total energy
        E = N / D
        E0 = E + mol.energy_nuc()

        #gradient
        Gpq = (Npq / D) - E * (Dpq / D)
        Gpq = Gpq - Gpq.conj().T
        Gvec = Gpq.flatten()
        Gnorm = lalg.norm(Gvec)

        timeEnergyGradient = time.time()

        #scipy optimizer
        paramsMat = - dt * Gpq[0:ne, ne:nso]
        params = paramsMat.flatten()
        params = complex_to_real(params)
        #params = np.zeros((2 * ne * (nso - ne), ), dtype = float)

        iterTol = tol
        if m == 0:
            iterTol = 1.e-6
        sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'SLSQP', jac = jac, tol = iterTol)
        #sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'SLSQP', tol = tol)
        #sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'SLSQP', jac = jac, tol = tol, options = {'maxiter' : 10})

        update = real_to_complex(sol.x)
        updateMat = update.reshape((ne, (nso - ne)))

        timeOptimizer = time.time()

        #build kappa matrix from params
        K = np.zeros(((nso, nso)), dtype = complex)
        K[0:ne, ne:nso] = updateMat
        K = K - K.conj().T

        #if optimizer fails use the gradient to update
        if sol.success == False:
            K = - dt * Gpq

        #orbital rotation
        U = lalg.expm(- K)

        #update orbitals
        orbs = orbs.dot(U)
        orbs = np.real(orbs) #when performing only Sz projection, we want real orbitals

        #update error
        error = abs(E - Eold)
        Eold = E

        #print
        if doPrint == True:
            print(f"-------------------------------- {m} --------------------------------")

            print("Projected values")
            print(f"  Denominator: {D}")
            print(f"  Numerator: {N}")
            print(f"  Electronic Energy: {E}")
            print(f"  Energy: {E0}")
            #print("Gradient")
            #print(np.real(Gpq))
            #print(np.imag(Gpq))
            print(f"  Gradient Norm: {Gnorm}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

            print("Scipy Optimizer")
            print(f"  message: " + sol.message)
            print(f"  fun: {sol.fun}")
            print(f"  jac: {lalg.norm(sol.jac)}")
            print(f"  nit: {sol.nit}")
            print(f"  Time for Optimizer: {timeOptimizer - timeEnergyGradient}")

            #print("Orbital Rotation")
            #print(U)
            #print("Orbitals")
            #print(orbs)

            print(f"Error: {error}")

        #check for convergence
        if (error < tol):
            break

    if doPrint == True:
        print(f"\nCalculation Complete")
        print(f"  Total Time: {time.time() - calcStart}")
        print(f"  Total Energy: {E0}")
    return E0, orbs


np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)

N = 30
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
#print("\n")

E0, mo = pGHF(mol, mf.mo_coeff)

#print("\n")
#print("Final Result")
#print("Energy")
#print(E0)
#print("Molecular Orbitals")
#print(mo)
