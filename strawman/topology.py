from typing import Dict, List, Optional, Tuple

import networkx
from networkx.algorithms.isomorphism import GraphMatcher
from openff.units import unit
from pydantic import Field, validator
from simtk.openmm import app

from strawman.molecule import Angstrom, Conformer, GraphMolecule, Molecule
from strawman.utilities.models import IndexRange, InnerModel, WrappedModel
from strawman.utilities.typing import ArrayQuantity


class MoleculeInstance(Molecule):
    """An instance of a particular molecule type in a topology."""

    @property
    def index(self) -> int:
        """The index of this molecules in the full topology"""
        return self._index

    @property
    def atom_indices(self) -> IndexRange:
        """The range of indices of the atoms in the full topology."""
        return self._atom_indices

    def __init__(self, molecule_type: Molecule, index: int, atom_indices: IndexRange):
        """

        Args:
            molecule_type: The type of molecule that this instance corresponds to.
            index: The index of this molecules in the full topology.
            atom_indices: The range of indices of the atoms in the full topology.
        """
        super().__init__()
        self._inner_data = molecule_type._inner_data

        self._index = index
        self._atom_indices = atom_indices


class Topology(WrappedModel):
    class Data(InnerModel):
        """The inner data model for a topology object."""

        molecule_types: Tuple[GraphMolecule.Data, ...] = Field(
            tuple(),
            description="A definition of the types of molecule contained in the "
            "topology.",
        )
        instances: Dict[str, int] = Field(
            {},
            description="The ranges of indices which correspond to a "
            "particular molecule type. These indices must be contiguous, "
            "non-overlapping and span between 0 and n_total_molecules.",
        )

        cell_vectors: Optional[ArrayQuantity[Angstrom]] = Field(
            None,
            description="The cell vectors [Å] of the box containing the system with "
            "shape=(3, 3).",
        )
        coordinates: Optional[Conformer] = Field(
            None,
            description="The XYZ coordinates [Å] of each atom contained within the "
            "topology stored in a numpy array with shape=(n_atoms, 3)",
        )

        @validator("instances")
        def validate_instances(cls, value: Dict[str, int]):

            # Sort the index ranges.
            sorted_ranges = sorted(value.keys(), key=lambda x: int(x.split("-")[0]))

            # Validate the range types and make sure the ranges are contiguous.
            sorted_value = {}
            previous_end: Optional[int] = None

            for index_range in sorted_ranges:

                molecule_type = value[index_range]

                assert isinstance(
                    index_range, str
                ), "each key must be an IndexRange object."
                assert isinstance(
                    molecule_type, int
                ), "each value must be the integer index of a molecule type."

                start, end = index_range.split("-")

                start = int(start)
                end = int(end)

                assert end > start

                if previous_end is None:
                    assert start == 0, "the instances index ranges do not start from 0."
                else:
                    assert (
                        start == previous_end
                    ), "the instances index ranges are not contiguous."

                sorted_value[index_range] = molecule_type
                previous_end = end

            return sorted_value

    @property
    def molecule_types(self) -> Tuple[GraphMolecule, ...]:
        return self._inner_data.molecule_types

    @property
    def n_molecules(self) -> int:

        if len(self._inner_data.instances) == 0:
            return 0

        return int([*self._inner_data.instances][-1].split("-")[1])

    @property
    def molecules(self):

        ranges: List[IndexRange] = [
            IndexRange.from_key(key) for key in self._inner_data.instances
        ]

        atom_end_index = 0
        range_index = 0

        for index in range(self.n_molecules):

            if index >= int(ranges[range_index].end):
                range_index += 1

            index_range = ranges[range_index]

            type_index = self._inner_data.instances[index_range.to_key()]
            molecule_type = self.molecule_types[type_index]

            n_atoms = len(molecule_type.atoms)

            atom_start_index = atom_end_index
            atom_end_index += n_atoms

            yield MoleculeInstance(
                molecule_type=Molecule.from_dict(
                    {
                        **molecule_type.dict(),
                        "conformer": None
                        if self._inner_data.coordinates is None
                        else self._inner_data.coordinates[
                            atom_start_index:atom_end_index, :
                        ],
                    }
                ),
                index=index,
                atom_indices=IndexRange(start=atom_start_index, end=atom_end_index),
            )

    @property
    def cell_vectors(self) -> Optional[unit.Quantity]:

        if self._inner_data.cell_vectors is None:
            return None

        return self._inner_data.cell_vectors * unit.angstrom

    @property
    def coordinates(self) -> Optional[unit.Quantity]:

        if self._inner_data.coordinates is None:
            return None

        return self._inner_data.coordinates * unit.angstrom

    def molecule_at_index(self, index: int) -> MoleculeInstance:

        found_type_index = -1

        atom_start_index = 0
        atom_end_index = -1

        # We can assume ``instances`` is sorted and contiguous thanks to the custom
        # validator.
        for index_range_key, type_index in self._inner_data.instances.items():

            index_range = IndexRange.from_key(index_range_key)

            n_atoms = len(self.molecule_types[type_index].atoms)

            if index < index_range.start or index >= index_range.end:

                atom_start_index += (index - index_range.start) * n_atoms
                atom_end_index = atom_start_index + n_atoms

                found_type_index = type_index
                break

            atom_start_index += (index_range.end - index_range.start) * n_atoms

        if found_type_index < 0:
            raise KeyError("Index out of range.")

        return MoleculeInstance(
            molecule_type=Molecule.from_dict(
                {
                    **self.molecule_types[found_type_index].dict(),
                    "conformer": self._inner_data.coordinates[
                        atom_start_index:atom_end_index, :
                    ],
                }
            ),
            index=index,
            atom_indices=IndexRange(start=atom_start_index, end=atom_end_index),
        )

    @classmethod
    def from_openmm(
        cls, openmm_topology: app.Topology, molecule_types: List[GraphMolecule]
    ) -> "Topology":

        # Convert all of the unique molecules to graphs
        type_graphs = []

        for molecule_type in molecule_types:

            molecule_graph = networkx.Graph()

            for i, atom in enumerate(molecule_type.atoms):
                molecule_graph.add_node(i, atomic_number=atom.atomic_number)

            for bond in molecule_type.bonds:
                molecule_graph.add_edge(
                    bond.index_a, bond.index_b, bond_order=bond.formal_order
                )

            type_graphs.append(molecule_graph)

        # Convert all openMM mols to graphs
        openmm_graph = networkx.Graph()

        for atom in openmm_topology.atoms():
            openmm_graph.add_node(atom.index, atomic_number=atom.element.atomic_number)

        for bond in openmm_topology.bonds():
            openmm_graph.add_edge(
                bond.atom1.index, bond.atom2.index, bond_order=bond.order
            )

        # For each connected subgraph (molecule) in the topology, find its matching type
        type_indices = []

        for openmm_sub_graph in (
            openmm_graph.subgraph(c).copy()
            for c in networkx.connected_components(openmm_graph)
        ):

            matches = [
                i
                for i, type_graph in enumerate(type_graphs)
                if GraphMatcher(
                    openmm_sub_graph,
                    type_graph,
                    node_match=lambda x, y: x["atomic_number"] == y["atomic_number"],
                    edge_match=lambda x, y: x["bond_order"] == y["bond_order"],
                ).is_isomorphic()
            ]

            if len(matches) != 1:
                raise RuntimeError("Wrong matches.")

            type_indices.append(matches[0])

        # Split the indices into ranges
        range_start = 0
        range_end = 0

        current_type = type_indices[0]

        instances = {}

        while len(type_indices) > 0:

            type_index = type_indices.pop(0)

            if type_index == current_type:

                range_end += 1
                continue

            instances[f"{range_start}-{range_end}"] = current_type

            range_start = range_end
            current_type = type_index

        if range_end > 0:

            range_end = range_end if range_end != range_start else range_end + 1
            instances[f"{range_start}-{range_end}"] = current_type

        return_value = cls()
        return_value._inner_data = cls.Data(
            molecule_types=tuple(
                GraphMolecule.Data.parse_obj(molecule_type.to_dict())
                for molecule_type in molecule_types
            ),
            instances=instances,
        )

        return return_value
