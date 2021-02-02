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
