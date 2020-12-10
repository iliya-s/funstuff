from pyscf import gto, scf, tools
import numpy as np
import scipy.linalg as lalg
import scipy.optimize as opt
import time
import sys

def contractEri(v2, dm1):
    nso = dm1.shape[0]
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

    return G1

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

def calcEnergy(S, H1, v2, Psi, Wg, Rg):
    Bra = Psi.copy()
    Ket = Psi.copy()

    D = 0.0
    N = 0.0
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket)

        #overlap quantities
        O = Bra.conj().T.dot(S).dot(Ketg)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #density matrix
        dm1 = Ketg.dot(invO).dot(Bra.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, dm1)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    E = N / D
    return np.real(E)

def calcGradient(S, H1, v2, Psi, Wg, Rg):
    Bra = Psi.copy()
    Ket = Psi.copy()

    D = 0.0
    N = 0.0
    Dab = np.zeros(Psi.shape, dtype = complex)
    Nab = np.zeros(Psi.shape, dtype = complex)
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket)

        #overlap quantities
        O = Bra.conj().T.dot(S).dot(Ketg)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to a,b orbital coefficient
        detOab  = detO * np.einsum('ai,i,ij,jb->ab', S, np.diag(Rg[i]), Psi, invO, dtype = complex)
        #detOab += detO * np.einsum('bi,ji,ja,a->ab', invO, Psi, S, np.diag(Rg[i]), dtype = complex)
        detOab = detOab + detOab

        #density matrix
        dm1 = Ketg.dot(invO).dot(Bra.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, dm1)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex)

        #drivative of H with respect to a,b orbital coefficient
        M = H1 + G1
        A1ab = np.einsum('bi,pi,pa,a->ab', invO, Psi, M, np.diag(Rg[i]), dtype = complex)
        #A2ab = np.einsum('aq,q,qi,ib->ab', M, np.diag(Rg[i]), Psi, invO, dtype = complex)
        A2ab = A1ab

        #B1ab = - np.einsum('ai,i,ij,jk,pk,pq,q,ql,lb->ab', S, np.diag(Rg[i]), Psi, invO, Psi, M, np.diag(Rg[i]), Psi, invO, dtype = complex)
        C1ak = np.einsum('ai,i,ij,jk->ak', S, np.diag(Rg[i]), Psi, invO, dtype = complex)
        D1pb = np.einsum('pq,q,ql,lb->pb', M, np.diag(Rg[i]), Psi, invO, dtype = complex)
        B1ab = - np.einsum('ak,pk,pb->ab', C1ak, Psi, D1pb)

        #B2ab = - np.einsum('bi,pi,pq,q,qj,jk,lk,la,a->ab', invO, Psi, M, np.diag(Rg[i]), Psi, invO, Psi, S, np.diag(Rg[i]), dtype = complex)
        #C2bq = np.einsum('bi,pi,pq,q->bq', invO, Psi, M, np.diag(Rg[i]), dtype = complex)
        #D2ja = np.einsum('jk,lk,la,a->ja', invO, Psi, S, np.diag(Rg[i]), dtype = complex)
        #B2ab = - np.einsum('bq,qj,ja->ab', C2bq, Psi, D2ja, dtype = complex)
        B2ab = B1ab

        Hab = A1ab + A2ab + B1ab + B2ab

        #averages with symmetry weights
        D += Wg[i] * detO
        Dab += Wg[i] * detOab

        N += Wg[i] * detO * H
        Nab += Wg[i] * (detOab * H + detO * Hab)

    E = N / D
    J = Nab / D - E * (Dab / D)
    return np.real(J)

def calcGradientFiniteDifference(S, H1, v2, Psi, Wg, Rg):
    ds = 1.e-6
    J = np.zeros(Psi.shape, dtype = complex)

    E0 = calcEnergy(S, H1, v2, Psi, Wg, Rg)
    for a in range(Psi.shape[0]):
        for b in range(Psi.shape[1]):
            Ket = Psi.copy()
            Ket[a, b] += ds
            E = calcEnergy(S, H1, v2, Ket, Wg, Rg)

            J[a, b] = (E - E0) / ds

    return J

#the following four functions are helpers to use scipy's optimization library

def real_to_complex(z): #real vector of length 2n -> complex of length n
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z): #complex vector of length n -> real of length 2n
    return np.concatenate((np.real(z), np.imag(z)))

def fun(params, nelectron, S, H1, v2, Wg, Rg):
    #wavefunction
    nso = S.shape[0]
    Psi = params.reshape((nso, nelectron))
    return calcEnergy(S, H1, v2, Psi, Wg, Rg)

def jac(params, nelectron, S, H1, v2, Wg, Rg):
    #wavefunction
    nso = S.shape[0]
    Psi = params.reshape((nso, nelectron))
    J = calcGradient(S, H1, v2, Psi, Wg, Rg)
    return J.flatten()

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

    assert(mo_coeff.shape[0] == nso)
    assert(mo_coeff.shape[1] == nso)
    #Psi will always be the wavefunction
    Psi = mo_coeff[:, 0:ne]

    #slight amount of noise helps optimization
    Psi += np.random.randn(nso, ne) / 100

    #Sz symmetry projector
    nGrid = 10
    Wg, Rg = calcSzSymmetryProjector(nao, sz, nGrid)

    Eold = 100
    dt = 0.1
    tol = 1.e-8
    doPrint = True
    calcStart = time.time()
    for m in range(20):
        iterStart = time.time()

        #energy
        E = calcEnergy(S, H1, v2, Psi, Wg, Rg)

        #gradient
        J = calcGradient(S, H1, v2, Psi, Wg, Rg)
        #Jfd = calcGradientFiniteDifference(S, H1, v2, Psi, Wg, Rg)

        #total energy
        E0 = E + mol.energy_nuc()

        #gradient
        Jvec = J.flatten()
        Jnorm = lalg.norm(Jvec)

        timeEnergyGradient = time.time()

        #scipy optimizer
        params = Psi.flatten()
        #params = Psi.flatten() - dt * Jvec

        sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'SLSQP', jac = jac, tol = tol)
        #sol = opt.minimize(fun, params, args = (ne, S, H1, v2, Wg, Rg), method = 'L-BFGS-B', jac = jac, tol = tol)

        #update parameters
        Psi = sol.x.copy().reshape((nso, ne))
        #if sol.success == True:
        #    Psi = sol.x.reshape((nso, ne))
        #else:
        #    Psi = Psi - dt * J

        timeOptimizer = time.time()

        #update error
        error = abs(E - Eold)
        Eold = E

        #print
        if doPrint == True:
            print(f"-------------------------------- {m} --------------------------------")

            print("Projected values")
            print(f"  Electronic Energy: {E}")
            print(f"  Total Energy: {E0}")
            #print("Gradient")
            #print(np.real(Gpq))
            #print(np.imag(Gpq))
            print(f"  Gradient Norm: {Jnorm}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

            print("Scipy Optimizer")
            print(f"  message: {sol.message}")
            print(f"  fun: {sol.fun}")
            print(f"  jac: {lalg.norm(sol.jac)}")
            print(f"  nit: {sol.nit}")
            print(f"  Time for Optimizer: {timeOptimizer - timeEnergyGradient}")

            #print("Occupied Orbitals")
            #print(Psi)

            print(f"Error: {error}")

        #check for convergence
        if (error < tol):
            break

    if doPrint == True:
        print(f"\nCalculation Complete")
        print(f"  Total Time: {time.time() - calcStart}")
        print(f"  Total Energy: {E0}")
    return E0, Psi


np.set_printoptions(precision=6)
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
#print(fock)
#print("\n")

E0, mo = pGHF(mol, mf.mo_coeff)

#print("\n")
#print("Final Result")
#print("Energy")
#print(E0)
#print("Molecular Orbitals")
#print(mo)
