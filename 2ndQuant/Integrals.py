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

'''
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

I2 = mol.intor('int2e',aosym='s8')
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
'''
