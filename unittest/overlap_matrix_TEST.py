#!/usr/bin/python3

from pyscf import gto
import json
import itertools

mol = gto.Mole()

mol.atom = "H 0 0 0; H 0 0 0.566918"
mol.basis = "sto-3g"
mol.spin = None
mol.build()

print("atom  : ", mol.atom)
print("_atom : ", mol._atom)
print("x: ", mol._atom[0][1][0])
print("y: ", mol._atom[0][1][1])
print("z: ", mol._atom[0][1][2])
print("basis : ", mol.basis)
print("_basis: ", mol._basis)

overlap = mol.intor("int1e_ovlp")
print("overlap: \n", overlap)