from strawman.molecule import BondStereochemistry, Molecule


def test_from_smiles():

    molecule = Molecule.from_smiles("C/C=C/C")

    assert molecule.to_smiles() == "C/C=C/C"
    assert molecule.atom_stereochemistry == {}
    assert molecule.bond_stereochemistry == {1: BondStereochemistry.Trans}

    with open("molecule.json", "w") as file:
        file.write(molecule.to_json())

    with open("molecule.json", "r") as file:
        loaded_molecule = Molecule.from_json(file.read())

    print(loaded_molecule)
