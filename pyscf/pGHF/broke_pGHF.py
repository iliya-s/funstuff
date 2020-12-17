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
    Jaa = np.einsum('pqrs,sr->pq', v2, dm1aa, dtype = complex, optimize = True)

    Jbb = np.einsum('pqrs,sr->pq', v2, dm1bb, dtype = complex, optimize = True)

    Kaa = np.einsum('psrq,sr->pq', v2, dm1aa, dtype = complex, optimize = True)
    Kbb = np.einsum('psrq,sr->pq', v2, dm1bb, dtype = complex, optimize = True)

    Kab = np.einsum('psrq,sr->pq', v2, dm1ba, dtype = complex, optimize = True)
    Kba = np.einsum('psrq,sr->pq', v2, dm1ab, dtype = complex, optimize = True)

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

def calcEnergy(S, H1, v2, Psi, Wg, Rg):
    Bra = Psi.copy()
    Ket = Psi.copy()

    #denominator and numerator of energy
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
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    #for complex conjugation projection, calculate overlaps with the conjugate of the ket
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket.conj())

        #overlap quantities
        O = Bra.conj().T.dot(S).dot(Ketg)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #density matrix
        dm1 = Ketg.dot(invO).dot(Bra.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, dm1)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        N += Wg[i] * detO * H

    E = N / D
    return np.real(E)

#this returns the gradient with respect to a,b and a,b* parameters of the energy
def calcComplexGradient(S, H1, v2, Psi, Wg, Rg):
    Bra = Psi.copy()
    Ket = Psi.copy()

    #denominator and numerator of energy
    D = 0.0
    N = 0.0
    #gradient with respect to a,b orbital parameter of the denominator and numerator of energy
    Dab = np.zeros(Psi.shape, dtype = complex)
    Nab = np.zeros(Psi.shape, dtype = complex)
    #gradient with respect to a,b* orbital parameter of the denominator and numerator of energy
    Dab_bar = np.zeros(Psi.shape, dtype = complex)
    Nab_bar = np.zeros(Psi.shape, dtype = complex)

    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket)

        #overlap quantities
        O = Bra.conj().T.dot(S).dot(Ketg)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to a,b orbital parameter
        detOab = detO * np.einsum('bi,ji,ja,a->ab', invO, Psi.conj(), S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOab_bar = detO * np.einsum('ai,i,ij,jb->ab', S, np.diag(Rg[i]), Psi, invO, dtype = complex, optimize = True)

        #density matrix
        dm1 = Ketg.dot(invO).dot(Bra.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, dm1)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex, optimize = True)

        #drivative of H with respect to a,b orbital parameter
        M = H1 + G1

        Aab = np.einsum('bi,pi,pa,a->ab', invO, Psi.conj(), M, np.diag(Rg[i]), dtype = complex, optimize = True)
        Baj = np.einsum('jk,lk,la,a->aj', invO, Psi.conj(), S, np.diag(Rg[i]), dtype = complex, optimize = True)
        Hab = Aab - np.einsum('qb,qj,aj->ab', Aab, Psi, Baj, dtype = complex, optimize = True)

        Aab_bar = np.einsum('aq,q,qi,ib->ab', M, np.diag(Rg[i]), Psi, invO, dtype = complex, optimize = True)
        Bak_bar = np.einsum('ai,i,ij,jk->ak', S, np.diag(Rg[i]), Psi, invO, dtype = complex, optimize = True)
        Hab_bar = Aab_bar - np.einsum('ak,pk,pb->ab', Bak_bar, Psi.conj(), Aab_bar, dtype =  complex, optimize = True)

        #averages with symmetry weights
        D += Wg[i] * detO
        Dab += Wg[i] * detOab
        Dab_bar += Wg[i] * detOab_bar

        N += Wg[i] * detO * H
        Nab += Wg[i] * (detOab * H + detO * Hab)
        Nab_bar += Wg[i] * (detOab_bar * H + detO * Hab_bar)

    #for complex conjugation projection, calculate overlaps with the conjugate of the ket, note that the only nonzero derivatives are with respect to a,b* parameters
    for i in range(len(Wg)):
        Ketg = Rg[i].dot(Ket.conj())

        #overlap quantities
        O = Bra.conj().T.dot(S).dot(Ketg)
        invO = lalg.inv(O)
        detO = lalg.det(O)

        #derivative of detO with respect to a,b orbital parameter
        detOab_bar  = detO * np.einsum('ai,i,ij,jb->ab', S, np.diag(Rg[i]), Psi.conj(), invO, dtype = complex, optimize = True)
        #detOab_bar += detO * np.einsum('bi,ji,ja,a->ab', invO, Psi.conj(), S, np.diag(Rg[i]), dtype = complex, optimize = True)
        detOab_bar = detOab_bar + detOab_bar

        #density matrix
        dm1 = Ketg.dot(invO).dot(Bra.conj().T)

        #hamiltonian quantities
        G1 = contractEri(v2, dm1)
        F1 = H1 + 0.5 * G1
        H = np.einsum('pq,qp->', F1, dm1, dtype = complex, optimize = True)

        #drivative of H with respect to a,b orbital parameter
        M = H1 + G1

        A1ab_bar = np.einsum('bi,pi,pa,a->ab', invO, Psi.conj(), M, np.diag(Rg[i]), dtype = complex, optimize = True)
        #A2ab_bar = np.einsum('aq,q,qi,ib->ab', M, np.diag(Rg[i]), Psi.conj(), invO, dtype = complex, optimize = True)
        A2ab_bar = A1ab_bar

        C1ak_bar = np.einsum('ai,i,ij,jk->ak', S, np.diag(Rg[i]), Psi.conj(), invO, dtype = complex, optimize = True)
        #D1pb_bar = np.einsum('pq,q,ql,lb->pb', M, np.diag(Rg[i]), Psi.conj(), invO, dtype = complex, optimize = True)
        D1pb_bar = A1ab_bar
        B1ab_bar = - np.einsum('ak,pk,pb->ab', C1ak_bar, Psi.conj(), D1pb_bar, optimize = True)

        #C2bq_bar = np.einsum('bi,pi,pq,q->bq', invO, Psi.conj(), M, np.diag(Rg[i]), dtype = complex, optimize = True)
        #D2ja_bar = np.einsum('jk,lk,la,a->ja', invO, Psi.conj(), S, np.diag(Rg[i]), dtype = complex, optimize = True)
        #B2ab_bar = - np.einsum('bq,qj,ja->ab', C2bq_bar, Psi.conj(), D2ja_bar, dtype = complex, optimize = True)
        B2ab_bar = B1ab_bar

        Hab_bar = A1ab_bar + A2ab_bar + B1ab_bar + B2ab_bar

        #averages with symmetry weights
        D += Wg[i] * detO
        Dab_bar += Wg[i] * detOab_bar

        N += Wg[i] * detO * H
        Nab_bar += Wg[i] * (detOab_bar * H + detO * Hab_bar)

    E = N / D
    J = Nab / D - E * (Dab / D)
    J_bar = Nab_bar / D - E * (Dab_bar / D)
    return J, J_bar

#this will calculate the gradient with respect to the real and imaginary part of the orbital parameters, given the gradients with respect to a,b and a,b*
def calcComplexGradient(J, J_bar):
    Jr = J + J_bar
    Ji = 1j * (J - J_bar)
    return np.real(Jr), np.real(Ji)

#U here is a given molecular orbital rotation
def calcNonRedundantEnergy(S, H1, v2, nelectron, mo_coeff, U, Wg, Rg):
    nso = S.shape[0]
    orbs = mo_coeff.dot(U)
    Psi = orbs[:, 0:nelectron]
    return calcEnergy(S, H1, v2, Psi, Wg, Rg)

def calcNonRedundantGradient(S, H1, v2, nelectron, mo_coeff, U, Wg, Rg):
    nso = S.shape[0]
    orbs = mo_coeff.dot(U)
    Psi = orbs[:, 0:nelectron]

    Jmi, Jmi_bar = calcComplexGradient(S, H1, v2, Psi, Wg, Rg)

    Gqi = np.einsum('mi,mq->qi', Jmi, mo_coeff, dtype = complex, optimize = True)
    Gqi_bar = np.einsum('mi,mq->qi', Jmi_bar, mo_coeff.conj(), dtype = complex, optimize = True)
    return calcGradient(Gqi, Gqi_bar)

#finite difference gradients
def calcGradientFiniteDifference(S, H1, v2, Psi, Wg, Rg):
    ds = 1.e-6
    E0 = calcEnergy(S, H1, v2, Psi, Wg, Rg)

    Jr = np.zeros(Psi.shape, dtype = float)
    for a in range(Psi.shape[0]):
        for b in range(Psi.shape[1]):
            Ket = Psi.copy()
            Ket[a, b] += ds
            E = calcEnergy(S, H1, v2, Ket, Wg, Rg)

            Jr[a, b] = (E - E0) / ds

    Ji = np.zeros(Psi.shape, dtype = float)
    for a in range(Psi.shape[0]):
        for b in range(Psi.shape[1]):
            Ket = Psi.copy()
            Ket[a, b] += 1j * ds
            E = calcEnergy(S, H1, v2, Ket, Wg, Rg)

            Ji[a, b] = (E - E0) / ds

    return Jr, Ji

def calcNonRedundantGradientFiniteDifference(S, H1, v2, nelectron, mo_coeff, U, Wg, Rg):
    ds = 1.e-6
    nso = S.shape[0]

    orbs = mo_coeff.dot(U)
    Psi = orbs[:, 0:nelectron]
    E0 = calcEnergy(S, H1, v2, Psi, Wg, Rg)

    Jr = np.zeros((nso, nelectron), dtype = float)
    for a in range(Jr.shape[0]):
        for j in range(Jr.shape[1]):
            Ucopy = U.copy()
            Ucopy[a, j] += ds
            orbs = mo_coeff.dot(Ucopy)
            Psi = orbs[:, 0:nelectron]
            E = calcEnergy(S, H1, v2, Psi, Wg, Rg)

            Jr[a, j] = (E - E0) / ds

    Ji = np.zeros((nso, nelectron), dtype = float)
    for a in range(Ji.shape[0]):
        for j in range(Ji.shape[1]):
            Ucopy = U.copy()
            Ucopy[a, j] += 1j * ds
            orbs = mo_coeff.dot(Ucopy)
            Psi = orbs[:, 0:nelectron]
            E = calcEnergy(S, H1, v2, Psi, Wg, Rg)

            Ji[a, j] = (E - E0) / ds

    return Jr, Ji

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

#the following four functions are helpers to use scipy's optimization library

def real_to_complex(z): #real vector of length 2n -> complex of length n
    return z[:len(z)//2] + 1j * z[len(z)//2:]

def complex_to_real(z): #complex vector of length n -> real of length 2n
    return np.concatenate((np.real(z), np.imag(z)))

def fun(params, S, H1, v2, nelectron, mo_coeff, Wg, Rg):
    #build parameter matrix, Uai block
    paramVec = real_to_complex(params)
    nso = S.shape[0]
    paramMat = paramVec.reshape((nso - nelectron, nelectron))
    #Unitary transformation
    U = np.identity(nso, dtype = complex)
    U[nelectron:nso, 0:nelectron] = paramMat
    U[0:nelectron, nelectron:nso] = - paramMat.conj().T
    return calcNonRedundantEnergy(S, H1, v2, nelectron, mo_coeff, U, Wg, Rg)

def jac(params, S, H1, v2, nelectron, mo_coeff, Wg, Rg):
    #build parameter matrix, Uai block
    paramVec = real_to_complex(params)
    nso = S.shape[0]
    paramMat = paramVec.reshape((nso - nelectron, nelectron))
    #Unitary transformation
    U = np.identity(nso, dtype = complex)
    U[nelectron:nso, 0:nelectron] = paramMat
    U[0:nelectron, nelectron:nso] = - paramMat.conj().T
    #gradient
    Jr, Ji = calcNonRedundantGradient(S, H1, v2, nelectron, mo_coeff, U, Wg, Rg)
    Jr_vec = Jr[nelectron:nso, 0:nelectron].flatten()
    Ji_vec = Ji[nelectron:nso, 0:nelectron].flatten()
    return np.concatenate((Jr_vec, Ji_vec), axis = 0)

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
    orbs = mo_coeff.copy()

    #randomly mix occ-occ and virt-virt orbitals via a complex transformation, slight amount of noise helps optimization
    randOccMat = randOrthoMat(ne)
    randVirMat = randOrthoMat(nso - ne)
    randU = lalg.block_diag(randOccMat, randVirMat)
    orbs = orbs.dot(randU)

    #initial orbital rotation is the identity
    U0 = np.identity(nso, dtype = complex)

    #Sz symmetry projector
    nGrid = 8
    Wg, Rg = calcSzSymmetryProjector(nso, sz, nGrid)

    Eold = 100
    tol = 1.e-8
    dt = 0.1
    doPrint = True
    calcStart = time.time()
    for m in range(50):
        iterStart = time.time()

        #initial orbital rotation
        U = U0.copy()

        #energy
        E = calcNonRedundantEnergy(S, H1, v2, ne, orbs, U, Wg, Rg)

        #gradient
        Jr, Ji = calcNonRedundantGradient(S, H1, v2, ne, orbs, U, Wg, Rg)
        print(Jr)
        print(Ji)
        Jfdr, Jfdi = calcNonRedundantGradientFiniteDifference(S, H1, v2, ne, orbs, U, Wg, Rg)
        print(Jfdr)
        print(Jfdi)

        #total energy
        E0 = E + mol.energy_nuc()

        #gradient
        J = Jr + 1j * Ji
        J_vec = J[ne:nso, 0:ne].flatten()
        J_norm = lalg.norm(J_vec)

        timeEnergyGradient = time.time()

        #gradient descent
        updateMat = - dt * J[ne:nso, 0:ne]

        #scipy optimizer
        paramsMat = np.zeros((nso - ne, ne), dtype = complex)
        paramsVec = paramsMat.flatten()
        params = complex_to_real(paramsVec)

        #params = params - dt * Jvec
        print(fun(params, S, H1, v2, ne, orbs, Wg, Rg))
        print(jac(params, S, H1, v2, ne, orbs, Wg, Rg))

        #sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'SLSQP', jac = jac, tol = tol)
        #sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'SLSQP')
        #sol = opt.minimize(fun, params, args = (S, H1, v2, ne, orbs, Wg, Rg), method = 'L-BFGS-B', jac = jac, tol = tol)

        #updateVec = real_to_complex(sol.x)
        #updateMat = updateVec.resize((nso - ne, ne))

        #new rotation matrix
        U[ne:nso, 0:ne] = updateMat
        U[0:ne, ne:nso] = - updateMat.conj().T

        #update parameters
        orbs = orbs.dot(U)

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
            print(f"  Gradient Norm: {J_norm}")
            print(f"  Time for Energy and Gradient: {timeEnergyGradient - iterStart}")

            #print("Scipy Optimizer")
            #print(f"  message: {sol.message}")
            #print(f"  fun: {sol.fun}")
            #print(f"  jac: {lalg.norm(sol.jac)}")
            #print(f"  nit: {sol.nit}")
            #print(f"  Time for Optimizer: {timeOptimizer - timeEnergyGradient}")

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
    return E0, orbs


np.set_printoptions(precision=4)
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
