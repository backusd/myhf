#!/usr/bin/python3

from pyscf import gto
import json
import itertools

mol = gto.Mole()
#mol.atom = "H 0 0 0; F 1 0.5 0.75"
mol.atom = "H 0 1.43233673 -0.96104039; H 0 -1.43233673 -0.96104039; O 0 0 0.2402601"
mol.basis = "sto-3g"
mol.spin = None
mol.verbose = 4
mol.unit = 'bohr'
mol.build()

print("atom  : ", mol.atom)
print("_atom : ", mol._atom)
print("x: ", mol._atom[0][1][0])
print("y: ", mol._atom[0][1][1])
print("z: ", mol._atom[0][1][2])
print("basis : ", mol.basis)
print("_basis: ", mol._basis)

# int1e_ovlp - overlap matrix
# int1e_kin  - kinetic energy matrix
# int1e_nuc  - nuclear electron attraction matrix
mat = mol.intor("int1e_nuc")
print("Matrix: \n", mat)

for row in mat:
    for val in row:
        print(val, end=' ')
    print()