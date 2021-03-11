from typing import Dict, List, Optional, Tuple

from pydantic import Field

from strawman.molecule import Atom, Bond, Conformer, GraphMolecule, Molecule
from strawman.utilities.models import ImmutableModel, WrappedModel


class Biomonomer(WrappedModel):
    class Data(ImmutableModel):

        identifier: Optional[str] = Field(
            None,
            description="An optional identifier associated with this monomer. In the "
            "case of polymers this would be a residue name, or in the case of DNA the "
            "code associated with a particular base.",
        )

        atoms: Tuple[Atom] = Field(
            ..., description="The atoms associated with this monomer."
        )

    @property
    def identifier(self) -> str:
        return self._inner_data.identifier

    @property
    def atoms(self) -> Tuple[Atom]:
        return self._inner_data.atoms


class GraphBiopolymer(GraphMolecule):
    class Data(ImmutableModel):

        identifier: Optional[str] = Field(
            None,
            description="An optional identifier associated with this biopolymer. This "
            "may be, for example, the chain id associated with a polypeptide.",
        )

        monomers: Tuple[Biomonomer] = Field(
            tuple(), description="The monomers which form the larger biopolymer."
        )
        bonds: Tuple[Bond] = Field(
            tuple(),
            description="The bonds between the atoms in the different monomers.",
        )

    @property
    def identifier(self) -> str:
        inner_data: GraphBiopolymer.Data = self._inner_data
        return inner_data.identifier

    @property
    def atoms(self) -> Tuple[Atom]:
        inner_data: GraphBiopolymer.Data = self._inner_data
        return tuple(atom for monomer in inner_data.monomers for atom in monomer.atoms)

    @property
    def monomers(self) -> Tuple[Biomonomer]:
        return tuple(self._inner_data.monomers)


class Biopolymer(GraphBiopolymer, Molecule):
    class Data(GraphBiopolymer.Data):

        conformer: Optional[Conformer] = Field(
            None,
            description="The XYZ coordinates [Å] of this molecules conformer stored "
            "in a numpy array with shape=(n_atoms, 3)",
        )

    @classmethod
    def from_molecule(
        cls, molecule: Molecule, residue_patterns: Optional[Dict[str, str]] = None
    ) -> "Biopolymer":

        if residue_patterns is None:

            residue_patterns = {
                "Ala": "[#6](-[#8])(=[#8])-[#6](-[#6])-[#7]",
                "Val": "[#7]-[#6](-[#6](-[#6])(-[#6]))-[#C](=[#8])",
            }

    @classmethod
    def from_monomers(cls, sequence: List[str]):
        pass
