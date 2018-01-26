from __future__ import division, print_function, absolute_import
import unittest
from rdkit import Chem
from rdkommtools import utils
from simtk import unit
from simtk.openmm import app
import numpy as np

class ConversionTester(unittest.TestCase):
    def test_rdmol_to_openmmTop(self):
        protein_fn = "tests/data/T4-protein.pdb"
        mol = Chem.MolFromPDBFile(protein_fn)

        top, omm_pos = utils.rdmol_to_openmmTop(mol)

        # Assert Atom numbers
        self.assertEqual(top.getNumAtoms(), mol.GetNumAtoms())

        for (op_at, rd_at) in zip(top.atoms(), mol.GetAtoms()):
            # Assert atom indexes
            self.assertEqual(op_at.index, rd_at.GetIdx())

        rd_pos = mol.GetConformer().GetPositions()
        # Assert atom positions
        np.testing.assert_almost_equal(np.array(rd_pos), omm_pos.in_units_of(unit.angstrom)/unit.angstrom)
        # self.assertEqual(rd_pos, omm_pos.in_units_of(unit.angstrom)/unit.angstrom)

        # Assert bond order
        # bond type is not matching, so do not compare
        dic_bond_openmm = {}
        for bond in top.bonds():
            # OpenMM atoms
            at0_idx = bond[0].index
            at1_idx = bond[1].index
            if at0_idx < at1_idx:
                dic_bond_openmm[(at0_idx, at1_idx)] = bond.order
            else:
                dic_bond_openmm[(at1_idx, at0_idx)] = bond.order

        dic_bond_rd = {}
        for bond in mol.GetBonds():
            # OE atoms
            at0_idx = bond.GetBeginAtom().GetIdx()
            at1_idx = bond.GetEndAtom().GetIdx()
            if at0_idx < at1_idx:
                dic_bond_rd[(at0_idx, at1_idx)] = bond.GetBondTypeAsDouble()
            else:
                dic_bond_rd[(at1_idx, at0_idx)] = bond.GetBondTypeAsDouble()

        self.assertEqual(dic_bond_openmm, dic_bond_rd)

    def test_openmmTop_to_rdmol(self):
        protein_fn ='tests/data/T4-protein.pdb'

        pdb = app.PDBFile(protein_fn)

        rd_mol = utils.openmmTop_to_rdmol(pdb.topology, pdb.positions)

        # Assert
        self.assertEqual(pdb.topology.getNumAtoms(), rd_mol.GetNumAtoms())

        for (op_at, rd_at) in zip(pdb.topology.atoms(), rd_mol.GetAtoms()):
            self.assertEqual(op_at.index, rd_at.GetIdx())

        rd_pos = rd_mol.GetConformer().GetPositions()
        np.testing.assert_almost_equal(pdb.getPositions(asNumpy=True).in_units_of(unit.angstrom)/unit.angstrom,
                                       np.array(rd_pos), decimal=2)
