from pyscf import gto, scf, tools
import numpy as np
import sys

#build molecule
mol = gto.Mole()
mol.atom = 'C 0 0 0; O 0 0 2.132'
mol.basis = 'sto-3g'
#mol.basis = 'sto-6g'
mol.verbose = 4
mol.spin = 0
mol.charge = 0
mol.unit = 'bohr'
mol.build()

#write atomic orbital integrals

t = mol.intor('int1e_kin')
v = mol.intor('int1e_nuc')
I1 = t + v
#t = mol.intor('int1e_kin_sph')
#v = mol.intor('int1e_nuc_sph')
print("t matrix")
print(t)
print("v matrix")
print(v)

I2 = mol.intor('int2e', aosym='s8')
#I2 = mol.intor('int2e_sph',aosym='s8')

tools.fcidump.from_integrals('AOFCIDUMP', I1, I2, mol.nao, mol.nelectron, mol.energy_nuc(), mol.spin)

s = mol.intor('int1e_ovlp')
#s = mol.intor('int1e_ovlp_sph')
f = open('METRIC','w')
for i in range(mol.nao):
    for j in range(mol.nao):
        f.write("%16.10e\t"% s[i, j])
    f.write("\n")
f.close()

mf = scf.RHF(mol)
mf.verbose=4
mf.kernel()
print("orbitals")
print(mf.mo_coeff)
