import pytest
from openff.toolkit.topology import Molecule as OFFMolecule
from openff.toolkit.topology import Topology as OFFTopology

from strawman.molecule import Molecule
from strawman.topology import Topology


def test_from_smiles():

    water = OFFMolecule.from_smiles("O")
    methane = OFFMolecule.from_smiles("C")

    openff_topology = OFFTopology.from_molecules(
        molecules=[water, water, methane, methane, water]
    )

    openmm_topology = openff_topology.to_openmm()

    straw_topology = Topology.from_openmm(
        openmm_topology, [Molecule.from_smiles("O"), Molecule.from_smiles("C")]
    )

    with open("topology.json", "w") as file:
        file.write(straw_topology.to_json())

    straw_molecules = [*straw_topology.molecules]
    print(openmm_topology)


@pytest.mark.xfail
def test_from_molecules_uniqueness():
    ethanol = Molecule.from_smiles("CCO")
    methane = Molecule.from_smiles("C")

    top = Topology()

    assert top.n_molecules == 0
    assert top.n_molecule_types == 0

    top.add_molecule(ethanol)
    top.add_molecule(methane)

    assert top.n_molecules == 2
    assert top.n_molecule_types == 2

    top.add_molecule(ethanol)

    assert top.n_molecules == 3
    assert top.n_molecule_types == 2


@pytest.mark.xfail
def test_combine_topologies():
    """Test the combination of heterogenous topologies"""

    ligand = Molecule.from_smiles("C1COC(=O)O1")
    water = Molecule.from_smiles("O")

    water_top = Topology.from_molecules(100 * [water])

    solv_top = water_top + ligand.to_topology()

    assert solv_top.n_molecules == 101
    assert solv_top.n_residues == 0  # Since none were specified (?)

    # Assert atom indices in ligand differ in the original molecule and the topology


@pytest.mark.xfail
def test_openmm_roundtrip():
    """Test the round-trip conversion through an OpenMM Topology object"""

    ethanol = Molecule.from_smiles("CCO")
    benzene = Molecule.from_smiles("c1ccccc1")
    top = Topology.from_molecules(molecules=[ethanol, benzene, benzene])

    omm_top = Topology.to_openmm()

    # TODO: Inspect the OpenMM Topology

    top_copy = Topology.from_openmm(omm_top, [ethanol, benzene])

    assert top.n_molecules == top_copy.n_molecules
    assert top.n_molecule_types == top_copy.n_molecule_types

    for mol, mol_copy in zip(top.molecules, top_copy.molecules):
        assert mol == mol_copy
