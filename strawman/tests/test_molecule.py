import pytest

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


@pytest.mark.xfail
def test_molecule_equality():

    ethane = Molecule.from_smiles("CC")
    ethane_copy = Molecule.from_smiles("CC")
    ethene = Molecule.from_smiles("C=C")
    chloroethane = Molecule.from_smiles("CCCl")
    z_chloropropene = Molecule.from_smiles("Cl\C=C/C")
    e_chloropropene = Molecule.from_smiles("Cl/C=C/C")

    assert ethane == ethane_copy
    assert ethane != ethene  # Bond order differs
    assert ethane != chloroethane  # Atom element differs
    assert z_chloropropene != e_chloropropene  # Bond stereochemistry differs

    # TODO: Atom stereochemistry differs

    # Could have an API point that allows more inspection that __eq__() ?

    ethane.is_equal_to(ethane_copy)

    with pytest.raises(BondMismatchError, match="order"):
        ethane.is_equal_to(ethene)

    with pytest.raises(AtomMistmatchError, match="element"):
        ethane.is_equal_to(chloroethane)

    with pytest.raises(BondMismatchError, match="stereo"):
        z_chloropropene.is_equal_to(e_chloropropene)
