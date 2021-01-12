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

    #1
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

    #2
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

    #3
    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.T, S, np.diag(Rg[i]), Phi, dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi, invO, Phi.T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    #4
    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.T, S, np.diag(Rg[i]), Phi.conj(), dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi.conj(), invO, Phi.T, dtype = complex, optimize = True)

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

    #1
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

    #2
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

        N += Wg[i] * (detO * H)
        Nmn_bar += Wg[i] * (detOmn_bar * H + detO * Hmn_bar)

    #3
    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.T, S, np.diag(Rg[i]), Phi, dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to m,n orbital parameter
        detOmn = detO * np.einsum('na,ab,bm,m->mn', invO, Phi.T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOmn += detO * np.einsum('ma,a,ab,bn->mn', S, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi, invO, Phi.T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #drivative of H with respect to m,n orbital parameter
        M1 = H1 + G1

        Amn1 = np.einsum('na,ab,bm,m->mn', invO, Phi.T, M1, np.diag(Rg[i]), dtype = complex, optimize = True)
        Amn2 = np.einsum('ma,a,ab,bn->mn', M1, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)

        Bmn1 = np.einsum('na,ab,bc,c,cd,de,ef,fm,m->mn', invO, Phi.T, M1, np.diag(Rg[i]), Phi, invO, Phi.T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        Bmn2 = np.einsum('ma,a,ab,bc,cd,de,e,ef,fn->mn', S, np.diag(Rg[i]), Phi, invO, Phi.T, M1, np.diag(Rg[i]), Phi, invO, dtype = complex, optimize = True)

        Hmn = Amn1 + Amn2 - Bmn1 - Bmn2

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        Dmn += Wg[i] * detOmn

        N += Wg[i] * detO * H
        Nmn += Wg[i] * (detOmn * H + detO * Hmn)

    #4
    for i in range(len(Wg)):
        #overlap quantities
        O = np.einsum('ma,ab,b,bn->mn', Phi.T, S, np.diag(Rg[i]), Phi.conj(), dtype = complex, optimize = True)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to m,n orbital parameter
        detOmn_bar = detO * np.einsum('na,ab,bm,m->mn', invO, Phi.T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOmn = detO * np.einsum('ma,a,ab,bn->mn', S, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)

        #transition density matrix
        tdm = np.einsum('m,ma,ab,bn->mn', np.diag(Rg[i]), Phi.conj(), invO, Phi.T, dtype = complex, optimize = True)

        #hamiltonian quantities
        G1 = contractEri(v2, tdm)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, tdm, dtype = complex, optimize = True)

        #drivative of H with respect to m,n orbital parameter
        M1 = H1 + G1

        Amn_bar = np.einsum('na,ab,bm,m->mn', invO, Phi.T, M1, np.diag(Rg[i]), dtype = complex, optimize = True)
        Bmn_bar = np.einsum('na,ab,bc,c,cd,de,ef,fm,m->mn', invO, Phi.T, M1, np.diag(Rg[i]), Phi.conj(), invO, Phi.T, S, np.diag(Rg[i]), dtype = complex, optimize = True)
        Hmn_bar = Amn_bar - Bmn_bar

        Amn = np.einsum('ma,a,ab,bn->mn', M1, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)
        Bmn = np.einsum('ma,a,ab,bc,cd,de,e,ef,fn->mn', S, np.diag(Rg[i]), Phi.conj(), invO, Phi.T, M1, np.diag(Rg[i]), Phi.conj(), invO, dtype = complex, optimize = True)
        Hmn = Amn - Bmn

        #averages with symmetry weights
        D += Wg[i] * detO
        Dmn += Wg[i] * detOmn
        Dmn_bar += Wg[i] * detOmn_bar

        N += Wg[i] * detO * H
        Nmn += Wg[i] * (detOmn * H + detO * Hmn)
        Nmn_bar += Wg[i] * (detOmn_bar * H + detO * Hmn_bar)

    E = N / D
    J = Nmn / D - E * (Dmn / D)
    J_bar = Nmn_bar / D - E * (Dmn_bar / D)

    Jr = np.real(J + J_bar)
    Ji = np.real(1j * (J - J_bar))

    J_vec = np.concatenate((Jr.flatten(), Ji.flatten()), axis = 0)
    return np.real(E), J_vec

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

    J_vec = np.concatenate((Jr.flatten(), Ji.flatten()), axis = 0)
    return E0, J_vec


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
    #calc energy
    E = calcEnergy(S, H1, v2, Phi, Wg, Rg)
    return E

def jac(params, nelectron, S, H1, v2, Wg, Rg):
    paramVec = real_to_complex(params)
    #wavefunction
    nso = S.shape[0]
    Phi = paramVec.reshape((nso, nelectron))
    #calc energy and gradient
    E, G = calcEnergyAndGradient(S, H1, v2, Phi, Wg, Rg)
    return G


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
    nGrid = 7
    Wg, Rg = calcSzSymmetryProjector(nso, sz, nGrid)

    #wavefunction
    Phi = orbs[:, 0:ne]
    #noise helps optimization
    Phi += 1j * 0.01 * (np.random.rand(Phi.shape[0], Phi.shape[1]) - 0.5)

    Eold = 100
    Etol = 1.e-8
    doPrint = True
    calcStart = time.time()
    for m in range(50):
        iterStart = time.time()

        #energy
        #E = calcEnergy(S, H1, v2, Phi, Wg, Rg)
        #print(E)

        #gradient
        E, G = calcEnergyAndGradient(S, H1, v2, Phi, Wg, Rg)
        #Efd, Gfd = calcEnergyAndGradientFiniteDifference(S, H1, v2, Phi, Wg, Rg)
        #print("G")
        #print(E)
        #print(G.T)
        #print("G fin-diff")
        #print(Efd)
        #print(Gfd.T)

        #total energy
        E0 = E + mol.energy_nuc()

        timeEnergyGradient = time.time()

        #calculate update

        #scipy optimizer
        params = complex_to_real(Phi.flatten())
        sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'SLSQP', jac = jac, tol = Etol)
        #sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'L-BFGS-B', jac = jac, tol = tol)

        #update parameters
        Phi = real_to_complex(sol.x).reshape((nso, ne))

        timeOptimizer = time.time()

        #update error
        Eerror = abs(E - Eold)
        Eold = E

        #print
        if doPrint == True:
            print(f"-------------------------------- {m} --------------------------------")

            print("Projected values")
            print(f"  Electronic Energy: {E}")
            print(f"  Total Energy: {E0}")
            print(f"  Gradient Norm: {lalg.norm(G)}")
            print(f"  Energy Error: {Eerror}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

            print("Scipy Optimizer")
            print(f"  message: {sol.message}")
            print(f"  fun: {sol.fun}")
            print(f"  jac: {lalg.norm(sol.jac)}")
            print(f"  nit: {sol.nit}")
            print(f"  Time for Optimizer: {timeOptimizer - timeEnergyGradient}")

        #check for convergence
        if Eerror < Etol:
            break

    if doPrint == True:
        print(f"\nCalculation Complete")
        print(f"  Total Time: {time.time() - calcStart}")
        print(f"  Total Energy: {E0}")
    return E0, Phi

np.set_printoptions(precision=2)
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

E0, mo = pGHF(mol, mf.mo_coeff)

#print("\n")
#print("Final Result")
#print("Energy")
#print(E0)
#print("Molecular Orbitals")
#print(mo)
