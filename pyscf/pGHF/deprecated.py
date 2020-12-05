from pyscf import gto, scf, tools
import numpy as np
import scipy.linalg as lalg
import scipy.optimize as opt
import sys

def calcLinearMethodMatrices(S, H1, V2, nelectron, mo_coeff, Wg, Rg):
    #hf det
    hf = [x for x in range(nelectron)]
    Psi = mo_coeff[:, hf]

    #singel excitations
    excitations, indices = generateSingleExcitations(nelectron, S.shape[0])
    nExcitations = len(excitations)

    A = np.zeros((1 + nExcitations, 1 + nExcitations), dtype = complex)
    B = np.zeros((1 + nExcitations, 1 + nExcitations), dtype = complex)

    #reference with reference
    N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi, Psi, Wg, Rg)

    A[0, 0] = N
    B[0, 0] = D

    #reference with tangent space
    for n in range(nExcitations):
        Psi_n = mo_coeff[:, indices[n]]
        N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi_n, Psi, Wg, Rg)

        A[0, 1 + n] = N
        A[1 + n, 0] = N.conj()

        B[0, 1 + n] = D
        B[1 + n, 0] = D.conj()

    #tangent space with tangent space
    for n in range(nExcitations):
        Psi_n = mo_coeff[:, indices[n]]

        N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi_n, Psi_n, Wg, Rg)

        A[1 + n, 1 + n] = N
        A[1 + n, 1 + n] = N.conj()

        B[1 + n, 1 + n] = D
        B[1 + n, 1 + n] = D.conj()

        for m in range(n):
            Psi_m = mo_coeff[:, indices[m]]

            N, D = calcEnergyNumeratorDenominator(S, H1, V2, Psi_n, Psi_m, Wg, Rg)

            A[1 + n, 1 + m] = N
            A[1 + m, 1 + n] = N.conj()

            B[1 + n, 1 + m] = D
            B[1 + m, 1 + n] = D.conj()

    return A, B

def calcFockMatrix(H1, V2, dm1):
    F = np.zeros(H1.shape, dtype = complex)
    F += H1.copy()
    F += np.einsum('pqrs,sr->pq', V2, dm1, casting = 'safe')
    F -= np.einsum('pqrs,qr->ps', V2, dm1, casting = 'safe')
    return F

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

