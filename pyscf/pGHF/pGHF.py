from pyscf import gto, scf, tools
import numpy as np
import scipy.linalg as lalg
import scipy.optimize as opt
import sys

def calcOverlap(S, Bra, Ket):
    ovlpMat = Bra.conj().T.dot(S).dot(Ket)
    return lalg.det(ovlpMat)

def calcDensityMatrix(S, Bra, Ket):
    ovlpMat = Bra.conj().T.dot(S).dot(Ket)
    ovlpInv = lalg.inv(ovlpMat)
    return Ket.dot(ovlpInv).dot(Bra.conj().T)

def calcHamiltonianMatrixElement(H1, V2, dm1):
    Hij = np.einsum('pq,qp->', H1, dm1, casting = 'safe')
    Hij += 0.5 * np.einsum('pqrs,qp,sr->', V2, dm1, dm1, casting = 'safe')
    Hij -= 0.5 * np.einsum('psrq,qp,sr->', V2, dm1, dm1, casting = 'safe')
    return Hij

def calcSymmetryProjector(nao, sz, nGrid = 8):
    Wg = []
    Rg = []
    for p in range(nGrid):
        theta = 2 * np.pi * p / nGrid

        w = np.exp(- 1j * theta * sz) / nGrid
        Wg.append(w)

        r = lalg.block_diag(np.exp(1j * theta * 0.5) * np.identity(nao), np.exp(1j * theta * (- 0.5)) * np.identity(nao))
        Rg.append(r)
    return Wg, Rg

def calcEnergyNumeratorDenominator(S, H1, V2, Bra, Ket, Wg, Rg):
    #integrate
    N = 0.0
    D = 0.0
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket)

        #calculate quantities
        Og = calcOverlap(S, Bra, Ketg)

        dmg = calcDensityMatrix(S, Bra, Ketg)
        Eg = calcHamiltonianMatrixElement(H1, V2, dmg)
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

def calcGradientEnergyNumeratorDenominator(S, H1, V2, nelectron, mo_coeff, Wg, Rg):
    Npq = np.zeros(S.shape, dtype = complex)
    Dpq = np.zeros(S.shape, dtype = complex)

    #hf det
    Psi = mo_coeff[:, 0:nelectron]

    #single excitations
    excitations, indices = generateSingleExcitations(nelectron, S.shape[0])
    for n in range(len(excitations)):
        Psi_n = mo_coeff[:, indices[n]]
        N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi_n, Psi, Wg, Rg)

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

def fun(params, S, H1, V2, nelectron, mo_coeff, Wg, Rg):
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
    N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi, Psi, Wg, Rg)
    N = np.real(N)
    D = np.real(D)
    return N / D

def jac(params, S, H1, V2, nelectron, mo_coeff, Wg, Rg):
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
    N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi, Psi, Wg, Rg)
    N = np.real(N)
    D = np.real(D)

    #gradient terms
    Npq, Dpq = calcGradientEnergyNumeratorDenominator(S, H1, V2, nelectron, orbs, Wg, Rg)

    G = (Npq / D) - (N / D) * (Dpq / D)
    G = G[0:nelectron, nelectron:nso]
    G = G.flatten()
    return complex_to_real(G)

def pGHF(mol, mo_coeff = None):
    #basic molecule info and integrals
    m = float(mol.spin) / 2
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
    V2 = np.zeros((nso, ) * 4)
    for i in range(nao):
        for j in range(nao):
            for k in range(nao):
                for l in range(nao):
                    ia = i
                    ja = j
                    ka = k
                    la = l

                    ib = i + nao
                    jb = j + nao
                    kb = k + nao
                    lb = l + nao

                    aaaa = (ia, ja, ka, la)
                    bbbb = (ib, jb, kb, lb)
                    aabb = (ia, ja, kb, lb)
                    bbaa = (ib, jb, ka, la)

                    V2[aaaa] = v2[(i, j, k, l)]
                    V2[bbbb] = v2[(i, j, k, l)]
                    V2[aabb] = v2[(i, j, k, l)]
                    V2[bbaa] = v2[(i, j, k, l)]

    if mo_coeff is None:
        orbs = np.random.randn(nso, nso)
        orbs = modifiedGramSchmidt(orbs, S)
    else:
        assert(mo_coeff.shape[0] == nso)
        assert(mo_coeff.shape[1] == nso)
        orbs = mo_coeff

    #symmetry projector
    Wg, Rg = calcSymmetryProjector(nao, m, 10)

    Eold = 100
    dt = 1.0
    tol = 1.e-10
    doPrint = True
    for m in range(50):
        #current wavefunction
        Psi = orbs[:, 0:ne]

        #energy terms
        N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi, Psi, Wg, Rg)
        N = np.real(N)
        D = np.real(D)

        #gradient terms
        Npq, Dpq = calcGradientEnergyNumeratorDenominator(S, H1, V2, ne, orbs, Wg, Rg)

        #electronic and total energy
        E = N / D
        E0 = E + mol.energy_nuc()

        #gradient
        G = (Npq / D) - E * (Dpq / D)
        G = G - G.conj().T
        Gvec = G.flatten()
        Gnorm = lalg.norm(Gvec)

        #scipy optimizer
        params = np.zeros((2 * ne * (nso - ne), ), dtype = float)
        sol = opt.minimize(fun, params, args = (S, H1, V2, ne, orbs, Wg, Rg), method = 'L-BFGS-B', jac = jac)
        update = real_to_complex(sol.x)
        paramMat = update.reshape((ne, (nso - ne)))
        K = np.zeros(((nso, nso)), dtype = complex)
        K[0:ne, ne:nso] = paramMat
        K = K - K.conj().T

        #orbital rotation
        #U = lalg.expm(dt * G)
        U = lalg.expm(- dt * K)
        #U = np.real(U) #when performing only Sz projection, we want real orbitals

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

            print(f"Denominator: {D}")
            print(f"Numerator: {N}")
            print(f"Electronic Energy: {E}")
            print(f"Energy: {E0}")

            print("Gradient")
            print(np.real(G))
            print(np.imag(G))

            print(f"Norm: {Gnorm}")

            print("\n scipy optimizer")
            print(sol)
            print("\n")

            print("Orbital Rotation")
            print(U)

            print("Orbitals")
            print(orbs)

            print(f"Error: {error}")

        #check for convergence
        if (error < tol):
            break

    return E, orbs


np.set_printoptions(precision=6)
np.set_printoptions(suppress=True)

mol = gto.Mole()
#mol.atom = 'H 0 0 0; H 0 0 1.4; H 0 0 2.8; H 0 0 4.2'
mol.atom = 'H 0 0 0; H 0 0 1.4; H 0 0 2.8;'
#mol.atom = 'H 0 0 0; H 0 0 1.4'
#mol.atom = 'H 0 0 0'
#mol.atom = 'Li 0 0 0'
#mol.basis = 'sto-3g'
#mol.basis = '631g'
#mol.basis = 'ccpvdz'
mol.verbose = 4
mol.spin = 1
mol.unit = 'bohr'
mol.build()

mf = scf.GHF(mol)
mf.kernel()
S = mf.get_ovlp()
occidx = mf.mo_occ > 0
occOrb = mf.mo_coeff[:, occidx]
print(occOrb)
fock = mf.get_hcore() + mf.get_veff()
print(fock)
print("\n")

E0, mo = pGHF(mol, mf.mo_coeff)
