from enum import Enum
from typing import Dict, List, Literal, Optional, Tuple

import numpy
from openeye import oechem, oeomega
from openff.units import unit
from pydantic import Field, conint

from strawman.utilities.elements import (
    AtomicNumber,
    AtomicSymbol,
    atomic_number_to_symbol,
)
from strawman.utilities.models import ImmutableModel, InnerModel, WrappedModel
from strawman.utilities.typing import ArrayQuantity

AtomIndex = conint(ge=0)
BondOrder = conint(ge=1, le=3)

Angstrom = Literal["angstrom"]

Conformer = ArrayQuantity[Angstrom]


class AtomStereochemistry(Enum):

    Undefined = "Undefined"

    Left = "Left"
    Right = "Right"

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.value


class BondStereochemistry(Enum):

    Undefined = "Undefined"

    Cis = "Cis"
    Trans = "Trans"

    def __str__(self):
        return self.value

    def __repr__(self):
        return self.value


OE_ATOM_STEREO_MAP = {
    oechem.OEAtomStereo_Undefined: AtomStereochemistry.Undefined,
    oechem.OEAtomStereo_Left: AtomStereochemistry.Left,
    oechem.OEAtomStereo_Right: AtomStereochemistry.Right,
}
OE_BOND_STEREO_MAP = {
    oechem.OEBondStereo_Undefined: BondStereochemistry.Undefined,
    oechem.OEBondStereo_Cis: BondStereochemistry.Cis,
    oechem.OEBondStereo_Trans: BondStereochemistry.Trans,
}


class Atom(ImmutableModel):

    # The core model fields which will be present in the serialized format.
    atomic_number: AtomicNumber = Field(
        ..., description="The atomic number of the atom."
    )
    formal_charge: int = Field(..., description="The formal charge on the atom.")

    # Derivative properties exposed via getters.
    @property
    def symbol(self) -> AtomicSymbol:
        return atomic_number_to_symbol(self.atomic_number)


class Bond(ImmutableModel):

    # The core model fields which will be present in the serialized format.
    index_a: AtomIndex = Field(..., description="The index of the bonds first atom.")
    index_b: AtomIndex = Field(..., description="The index of the bonds second atom.")

    formal_order: BondOrder = Field(
        ...,
        description="The formal bond order. This should be equal to 1 (single), "
        "2 (double), or 3 (triple) depending on the type of bond.",
    )
    fractional_order: Optional[float] = Field(
        None, description="The fractional bond order likely derived from a QM source."
    )


class GraphMolecule(WrappedModel):
    """A representation of a molecule which does not contain coordinates."""

    class Data(InnerModel):
        """The inner data model for a molecule object."""

        atoms: Tuple[Atom, ...] = Field(
            tuple(), description="The atoms which form the molecule."
        )
        bonds: Tuple[Bond, ...] = Field(
            tuple(), description="The bonds between the atoms."
        )

    @property
    def atoms(self) -> Tuple[Atom, ...]:
        return self._inner_data.atoms

    @property
    def bonds(self) -> Tuple[Bond, ...]:
        return self._inner_data.bonds

    def to_openeye(self) -> oechem.OEMol:

        oe_molecule = oechem.OEMol()

        # Add atoms and bonds.
        oe_atoms = []

        for atom in self.atoms:

            oe_atom = oe_molecule.NewAtom(atom.atomic_number)
            oe_atom.SetFormalCharge(atom.formal_charge)

            oe_atoms.append(oe_atom)

        for bond in self.bonds:

            oe_atom_a = oe_atoms[bond.index_a]
            oe_atom_b = oe_atoms[bond.index_b]

            oe_bond = oe_molecule.NewBond(oe_atom_a, oe_atom_b)
            oe_bond.SetOrder(bond.formal_order)

        oechem.OEFindRingAtomsAndBonds(oe_molecule)
        return oe_molecule

    def to_smiles(self) -> str:
        return oechem.OEMolToSmiles(self.to_openeye())


class Molecule(GraphMolecule):
    """A representation of a molecule which may contain coordinates."""

    class Data(GraphMolecule.Data):
        """The inner data model for a molecule object."""

        conformer: Optional[Conformer] = Field(
            None,
            description="The XYZ coordinates [Å] of this molecules conformer stored "
            "in a numpy array with shape=(n_atoms, 3)",
        )

    @property
    def conformer(self) -> Optional[unit.Quantity]:
        return self._inner_data.conformer * unit.angstrom

    @conformer.setter
    def conformer(self, value: Optional[Conformer]):
        self._inner_data.conformer = value

    @property
    def atom_stereochemistry(self) -> Dict[int, AtomStereochemistry]:
        return self._atom_stereochemistry()

    @property
    def bond_stereochemistry(self) -> Dict[int, BondStereochemistry]:
        return self._bond_stereochemistry()

    def _atom_stereochemistry(self) -> Dict[int, AtomStereochemistry]:
        """Adds a conformer to the molecule."""

        oe_molecule = self.to_openeye()

        atom_stereochemistry = {
            oe_atom.GetIdx(): oe_atom.GetStereo(
                [n for n in oe_atom.GetAtoms()], oechem.OEAtomStereo_Tetra
            )
            for oe_atom in oe_molecule.GetAtoms()
            if oe_atom.IsChiral()
        }

        return {
            index: OE_ATOM_STEREO_MAP[oe_stereo]
            for index, oe_stereo in atom_stereochemistry.items()
        }

    def _bond_stereochemistry(self) -> Dict[int, BondStereochemistry]:
        """Adds a conformer to the molecule."""

        oe_molecule = self.to_openeye()

        bond_stereochemistry = {
            oe_bond.GetIdx(): oe_bond.GetStereo(
                [
                    [n for n in oe_bond.GetBgn().GetAtoms() if n != oe_bond.GetEnd()][
                        0
                    ],
                    [n for n in oe_bond.GetEnd().GetAtoms() if n != oe_bond.GetBgn()][
                        0
                    ],
                ],
                oechem.OEBondStereo_CisTrans,
            )
            for oe_bond in oe_molecule.GetBonds()
            if oe_bond.IsChiral()
        }

        return {
            index: OE_BOND_STEREO_MAP[oe_stereo]
            for index, oe_stereo in bond_stereochemistry.items()
        }

    def to_openeye(self) -> oechem.OEMol:

        oe_molecule = super(Molecule, self).to_openeye()

        # Add a conformer if one is present and use it to perceive the
        # stereochemistry if possible.
        if self.conformer is not None:

            oe_molecule.DeleteConfs()
            oe_coords = oechem.OEFloatArray(self.conformer.magnitude.flatten())
            oe_molecule.NewConf(oe_coords)

            oechem.OEPerceiveChiral(oe_molecule)
            oechem.OE3DToInternalStereo(oe_molecule)

        oechem.OEFindRingAtomsAndBonds(oe_molecule)
        return oe_molecule

    @classmethod
    def from_openeye(cls, oe_molecule) -> List["Molecule"]:

        atoms: Dict[int, Atom] = {
            oe_atom.GetIdx(): Atom(
                atomic_number=oe_atom.GetAtomicNum(),
                formal_charge=oe_atom.GetFormalCharge(),
            )
            for oe_atom in oe_molecule.GetAtoms()
        }

        conformers = []

        for oe_conformer in oe_molecule.GetConfs():

            conformer = numpy.zeros((oe_molecule.NumAtoms(), 3))

            for atom_index, coordinates in oe_conformer.GetCoords().items():
                conformer[atom_index, :] = coordinates

            if oe_molecule.NumAtoms() > 1 and numpy.allclose(conformer, 0.0):
                continue

            conformers.append(conformer)

        if len(conformers) == 0:
            conformers.append(None)

        return_values = [
            cls.from_dict(
                Molecule.Data(
                    atoms=tuple(atoms[index] for index in range(len(atoms))),
                    bonds=tuple(
                        Bond(
                            index_a=oe_bond.GetBgnIdx(),
                            index_b=oe_bond.GetEndIdx(),
                            formal_order=oe_bond.GetOrder(),
                        )
                        for oe_bond in oe_molecule.GetBonds()
                    ),
                    conformer=conformer,
                ).dict()
            )
            for conformer in conformers
        ]

        return return_values

    @classmethod
    def from_smiles(cls, smiles: str) -> "Molecule":

        oe_molecule = oechem.OEMol()

        oechem.OESmilesToMol(oe_molecule, smiles)
        oechem.OEAddExplicitHydrogens(oe_molecule)

        # First create the base molecule object. This will ignore stereochemistry
        # in the smiles pattern.
        return_value = cls.from_openeye(oe_molecule)[0]

        # Next try to figure out if the molecule has stereochemistry defined
        # for all stereocenters.
        has_defined_stereochemistry = all(
            [
                *[
                    oe_atom.HasStereoSpecified()
                    for oe_atom in oe_molecule.GetAtoms()
                    if oe_atom.IsChiral()
                ],
                *[
                    oe_bond.HasStereoSpecified()
                    for oe_bond in oe_molecule.GetBonds()
                    if oe_bond.IsChiral()
                ],
            ]
        )

        # If the SMILES pattern did explicitly define stereochemistry try to generate
        # a conformer to encode this.
        if has_defined_stereochemistry:

            omega = oeomega.OEOmega()
            omega.SetMaxConfs(1)
            omega.SetIncludeInput(False)
            omega.SetEnergyWindow(15.0)
            omega.SetStrictStereo(True)
            omega.SetSampleHydrogens(True)
            omega.SetCanonOrder(False)
            omega(oe_molecule)

            for oe_conformer in oe_molecule.GetConfs():

                conformer = numpy.zeros((oe_molecule.NumAtoms(), 3))

                for atom_index, coordinates in oe_conformer.GetCoords().items():
                    conformer[atom_index, :] = coordinates

                return_value.conformer = conformer * unit.angstrom

        return return_value
