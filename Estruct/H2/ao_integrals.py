from pyscf import gto, scf, tools
import numpy as np
import sys

mol = gto.Mole()
#mol.atom = 'H 0 0 0; He 0 0 1.4632'
mol.atom = 'H 0 0 0; H 0 0 1.4'
#mol.basis = 'sto-6g'
mol.basis = 'sto-3g'
#mol.basis = 'ccpvdz'
mol.verbose = 4
#mol.charge = 1
mol.unit = 'bohr'
mol.build()

#s = mol.intor('int1e_ovlp')
#t = mol.intor('int1e_kin')
#v_nuc = mol.intor('int1e_nuc')
s = mol.intor('int1e_ovlp_sph')
t = mol.intor('int1e_kin_sph')
v_nuc = mol.intor('int1e_nuc_sph')
print(t)
print(v_nuc)

I1 = t + v_nuc
#I2 = mol.intor('int2e',aosym='s8')
I2 = mol.intor('int2e_sph',aosym='s8')

tools.fcidump.from_integrals('fcidump', I1, I2, mol.nao, mol.nelectron)

f = open('metric','w')
tools.dump_mat.dump_rec(f, s, ncol = mol.nao, digits = 14)
#tools.dump_mat.dump_rec(f, s, ncol = mol.nao, digits = 14, label = mol.ao_labels(), label2 = mol.ao_labels())
f.close()

f = open('e_nuc','w')
e_nuc = mol.energy_nuc()
f.write(str(e_nuc))
f.close()
