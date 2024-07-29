#!/usr/bin/python3

from pyscf import gto
import json
import itertools

atoms = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
        "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca"]

distances = [ 0.5 ]

basiss = [ "sto-3g", "sto-6g" ]

mol = gto.Mole()
results = []

for i in range(len(atoms)):
    for j in range(i, len(atoms)):
        atom1 = atoms[i]
        atom2 = atoms[j]

        for distance in distances:
            for basis in basiss:
                print("============================================")
                mol.atom = atom1 + " 0 0 0; " + atom2 + " " + str(distance) + " 0 0"
                mol.basis = basis
                mol.spin = None
                mol.build()

                print("atom  : ", mol.atom)
                print("_atom : ", mol._atom)
                print("basis : ", mol.basis)
                print("_basis: ", mol._basis)

                overlap = mol.intor("int1e_ovlp")
                print("overlap: ", overlap)

                for row in overlap:
                    for val in row:
                        print(val, end=' ')
                    print()

                result = {}
                result["atoms"] = [ atom1, atom2 ]
                result["positions"] = [mol._atom[0][1][0], mol._atom[0][1][1], mol._atom[0][1][2], mol._atom[1][1][0], mol._atom[1][1][1], mol._atom[1][1][2]]
                result["basis"] = basis
                result["overlap_rows"] = len(overlap)
                result["overlap_cols"] = len(overlap[0])
                vals = []
                for row in overlap:
                    for val in row:
                        vals.append(val.item())
                result["overlap"] = vals

                results.append(result)

resultsDict = { "results": results }

print(resultsDict)

json_object = json.dumps(resultsDict, indent=4)

# Writing results to file
with open("data/overlap-matrix-expected-values.json", "w") as outfile:
    outfile.write(json_object)

#mol.atom = '''O 0 0 0; H  0 1 0; H 0 0 1'''
#mol.atom = '''H  0 0 0; H 0.566918 0 0'''
#mol.basis = 'sto-3g'
#mol.build()

#print("Molecule:")

#d = vars(mol)
#for key in d:
#    print(key, ": ", d[key])

#print("atom  : ", mol.atom)
#print("_atom : ", mol._atom)
#print("basis : ", mol.basis)
#print("_basis: ", mol._basis)
#
#overlap = mol.intor("int1e_ovlp")
#print("overlap: ", overlap)
#
#for row in overlap:
#    for i in row:
#        print(i, end=' ')
#    print()
