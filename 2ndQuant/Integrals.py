from pyscf import gto, scf, tools
import numpy as np
import sys

#build molecule
#if H20
angle = 104.52 #Szabo
x = - np.sin((angle - 90) * np.pi / 180.0) * 1.809
y = np.cos((angle - 90) * np.pi / 180.0) * 1.809

mol = gto.Mole()

#mol.atom = 'H 0 0 0; H 0 0 1.4'
#mol.atom = 'H 0 0 0; He 0 0 1.4632'
mol.atom = 'Be 0 0 0'
#mol.atom = 'Mg 0 0 0'
#mol.atom = 'C 0 0 0'
#mol.atom = 'H %f %f 0; H 1.809 0 0;O 0 0 0' % (x, y)
#mol.atom = 'Ne 0 0 0'

mol.basis = 'sto-6g'
#mol.basis = '6-31g'
#mol.basis = 'ccpvdz'

mol.verbose = 5
mol.spin = 0
mol.charge = 0
mol.unit = 'bohr'
mol.build()

mf = scf.RHF(mol)
mf.kernel()
tools.fcidump.from_scf(mf, 'FCIDUMP')
