from pyscf import gto, scf, tools
import numpy as np
import sys

mol = gto.Mole()
mol.atom = 'H 0 0 0; H 0 0 1.4'
#mol.basis = 'sto-3g'
mol.basis = 'ccpvdz'
mol.verbose = 4
mol.spin = 0
mol.unit = 'bohr'
mol.build()

s = mol.intor('int1e_ovlp')
t = mol.intor('int1e_kin')
v_nuc = mol.intor('int1e_nuc')

I1 = t + v_nuc
I2 = mol.intor('int2e',aosym='s8')

tools.fcidump.from_integrals('fcidump', I1, I2, mol.nao, mol.nelectron)

f = open('s_matrix','w')
tools.dump_mat.dump_rec(f, s, ncol = mol.nao, digits = 14)
f.close()
