from pyscf import gto, scf, tools
import numpy as np
import sys

UHF = False

#build molecule
mol = gto.Mole()
mol.atom = "H 0 0 0; H 0 0 1.4"
#mol.atom = 'O 0 0 0; C 0 0 2.132'
#mol.basis = 'sto-3g'
mol.basis = 'sto-6g'
#mol.basis = 'ccpvdz'
mol.verbose = 4
mol.spin = 0
mol.unit = 'bohr'
mol.build()

#run scf
if UHF:
    mf = scf.UHF(mol)
else:
    mf = scf.RHF(mol)
mf.scf()

#write fcidump
print(mf.mo_coeff)
tools.fcidump.from_scf(mf, "fcidump")
