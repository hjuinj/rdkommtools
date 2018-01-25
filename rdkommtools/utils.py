from __future__ import division, absolute_imports, print_function
from rdkit import Chem, Geometry
from simtk.openmm import app, Vec3
from simtk import unit
import itertools
import pyparsing as pyp
from itertools import groupby

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

def rdmol_to_openmmTop(mol, confId = 0):
        """
    This function converts an rdmol to an openmm topology
    The rdmol coordinates are assumed to be in Angstrom unit
    Parameters:
    -----------
    mol: rdmol molecule
        The molecule to convert
    confId: int
        The id of the conformer from which coordinates will be taken from `mol`
    Return:
    -------
    topology : OpenMM Topology
        The generated OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the
        generated topology in Angstrom units
    """
    mol = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(mol))# hacky way to get all atom have PDB file relevant fields

    topology = app.Topology()
    rdk_atom_to_openmm = {}
    _chains = set([])

    atoms_grouped = groupby(mol.GetAtoms(), lambda atm :
    (atm.GetPDBResidueInfo().GetChainId(),
    atm.GetPDBResidueInfo().GetResidueNumber(),
    atm.GetPDBResidueInfo().GetResidueName())) #CESHI fails when not read from PDB
    for key, residue in atoms_grouped:
        chainId, resNum, resName = key
        if chainId not in _chains:
            _chains.add(chainId)
            openmm_chain = topology.addChain(chainId)

        openmm_res = topology.addResidue(resName, openmm_chain)
        for atm in residue:
            element = app.element.Element.getByAtomicNumber(atm.GetAtomicNum())
            openmm_at = topology.addAtom(atm.GetPDBResidueInfo().GetName().strip() , element, openmm_res)
            openmm_at.index = atm.GetIdx()
            rdk_atom_to_openmm[atm.GetIdx()] = openmm_at

    if topology.getNumAtoms() != mol.GetNumAtoms():
        raise ValueError("OpenMM topology and RDMol number of atoms mismatching: "
                      "OpenMM = {} vs RDMol  = {}".format(topology.getNumAtoms(), mol.GetNumAtoms()))


    # Count the number of bonds in the openmm topology
    omm_bond_count = 0

    def IsAmideBond(rdk_bond):
        # This supporting function checks if the passed bond is an amide bond or not.
        # Our definition of amide bond C-N between a Carbon and a Nitrogen atom is:
        #          O
        #          ║
        #  CA or O-C-N-
        #            |

        # The amide bond C-N is a single bond
        if str(rdk_bond.GetBondType()) != "SINGLE" :
            return False

        atomB, atomE = rdk_bond.GetBeginAtom(), rdk_bond.GetEndAtom()
        # The amide bond is made by Carbon and Nitrogen atoms
        if not (atomB.GetAtomicNum() == 6 and atomE.GetAtomicNum() == 7 or
        (atomB.GetAtomicNum() == 7 and atomE.GetAtomicNum() == 6)):
            return False

        # Select Carbon and Nitrogen atoms
        if atomB.GetAtomicNum() === 6 :
            C_atom = atomB
            N_atom = atomE
        else:
            C_atom = atomE
            N_atom = atomB

        # Carbon and Nitrogen atoms must have 3 neighbour atoms
        if not (C_atom.GetDegree() == 3 and N_atom.GetDegree() == 3): #CESHI
            return False

        double_bonds, single_bonds = 0, 0

        for bond in C_atom.GetBonds():
            # The C-O bond can be single or double.
            if (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 8 ) or (bond.GetBeginAtom().GetAtomicNum() == 8 and bond.GetEndAtom().GetAtomicNum() == 6 ):
                if str(bond.GetBondType()) == "DOUBLE":
                    double_bonds += 1
                if str(bond.GetBondType()) == "SINGLE":
                    single_bonds += 1
            # The CA-C bond is single
            if (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6 ):
                if str(bond.GetBondType()) == "SINGLE":
                    single_bonds += 1
        # Just one double and one single bonds are connected to C
        # In this case the bond is an amide bond
        if double_bonds == 1 and single_bonds == 1:
            return True
        else:
            return False

    # Creating bonds
    for bond in mol.GetBonds():
        omm_bond_count += 1

        # Set the bond type
        bond_order = bond.GetBondTypeAsDouble()
        if IsAmideBond(bond):
            omm_bond_type = "Amide"
        elif bond_order == 1.0:
            omm_bond_type = "Single"
        elif bond_order == 2.0:
            omm_bond_type = "Double"
        elif bond_order == 3.0:
            omm_bond_type = "Triple"
        elif bond_order == 1.5:
            omm_bond_type = "Aromatic"
        else:
            omm_bond_type = None

        topology.AddBond(
            rdk_atom_to_openmm[bond.GetBeginAtom().GetIdx()],
            rdk_atom_to_openmm[bond.GetEndAtom().GetIdx()],
            type = omm_bond_type,
            order = bond_order #CESHI the bond order calculated is a double, supposedly OpenMM takes an int value

    if omm_bond_count != mol.GetNumBonds():
        raise ValueError("OpenMM topology and RDMol number of bonds mismatching: "
                      "OpenMM = {} vs RDMol  = {}".format(omm_bond_count, mol.GetNumBonds()))
    coords = mol.GetConformer(confId).GetPositions()
    positions = [Vec3(v[0], v[1], v[2]) for v in coords] * unit.angstroms

    return topology, positions


def openmmTop_to_rdmol(topology, positions, verbose = False):
    """
    This function converts an OpenMM topology into a RDMol
    Parameters:
    -----------
    topology : OpenMM Topology
        The OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the
        topology
    Return:
    -------
    rdmol : RDMol
        The generated RDKit molecule
    """
    rdmol = Chem.RWMol()

    # Mapping dictionary between openmm atoms and rdk atoms
    openmm_atom_to_rdk_atom = {}

    # Python set used to identify atoms that are not in protein residues
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

    #TODO charge info is not transferred
    for chain in topology.chains():
        chainId = str(chain.id)
        for res in chain.residues():
            resName, resNum= res.name, int(res.idx)
            for openmm_at in res.atoms():
                rdatom = Chem.Atom(openmm_at.element._atomic_number)
                rdatom.GetPDBResidueInfo().SetName(openmm_at.name)
                rdatom.GetPDBResidueInfo().SetChainId(chainId)
                rdatom.GetPDBResidueInfo().SetResidueNumber(resNum)
                rdatom.GetPDBResidueInfo().SetResidueName(resName)

                if resName not in keep:
                    rdatom.SetIsHeteroAtom()

                rdmol.AddAtom(rdatom)
                openmm_atom_to_rdk_atom[openmm_at] = rdmol.GetNumAtoms() - 1

    if topology.getNumAtoms() != rdmol.GetNumAtoms():
        raise ValueError("OpenMM topology and RDMol number of atoms mismatching: "
                             "OpenMM = {} vs RDMol  = {}".format(topology.getNumAtoms(), rdmol.GetNumAtoms()))
    # Count the number of bonds in the openmm topology
    omm_bond_count = 0

    # Create the bonds
    _bondtypes = {0: Chem.BondType.UNSPECIFIED,
              1: Chem.BondType.SINGLE,
              1.5: Chem.BondType.AROMATIC,
              2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE,
              4: Chem.BondType.QUADRUPLE,
              5: Chem.BondType.QUINTUPLE,
              6: Chem.BondType.HEXTUPLE,
              7: Chem.BondType.ONEANDAHALF,}

    for omm_bond in topology.bonds():

        omm_bond_count += 1

        at0 = omm_bond[0]
        at1 = omm_bond[1]

        rd_atom0, rd_atom1 = openmm_atom_to_rdk_atom[at0], openmm_atom_to_rdk_atom[at1]

        if omm_bond.type == "Aromatic":
            #CESHI assumed by setting bond aromatic the two atoms are aromatic
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[1.5])
        elif omm_bond.type == "Single":
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[1])
        elif omm_bond.type == "Double":
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[2])
        elif omm_bond.type == "Triple":
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[3])
        elif omm_bond.type == "Amide":
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[int(omm_bond.order)])
        else:
            rdmol.AddBond(rd_atom0, rd_atom1, _bondtypes[0])

    if topology.getNumAtoms() != rdmol.GetNumAtoms():
        raise ValueError("OpenMM topology and RDMol number of bonds mismatching: "
                             "OpenMM = {} vs RDMol  = {}".format(omm_bond_count, rdmol.GetNumBonds()))

    pos = positions.in_units_of(unit.angstrom) / unit.angstrom
    conformer = Chem.Conformer()

    for idx,coord in enumerate(positions):
        x,y,z = [i._value for i in coord]
        conformer.SetAtomPosition(idx, Geometry.Point3D(x,y,z))

    rdmol.AddConformer(conformer)

    rdmol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(rdmol)

    return rdmol.GetMol()




def delete_shell(core_mol, del_mol, cut_off, in_out='in'):
    """
    This function deletes molecules present in the passed argument
    del_mol that are far (in_out=out) or close (in_out=in) than the
    selected cutoff distance (in A) from the passed molecules core_mol
    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecules
    del_mol: OEMol molecule
        The molecules to be deleted if their distances from the core_mol
        molecules are greater or closer that the selected cutoff distance
    cut_off: python float number
        The threshold distance in A used to mark atom for deletion
    in_out: python string
        A flag used to select if delete molecules far or close than
        the cutoff distance from the core_mol
    Return:
    -------
    reset_del: copy of del_mol where atoms have been deleted with
        reset atom indexes
    """



def check_shell(core_mol, check_mol, cutoff):
    """
    This function checks if at least one atomic distance from the passed
    check_mol molecule to the core_mol molecule is less than the selected
    cutoff distance in A.
    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecule
    check_mol: OEMol molecule
        The molecule to be checked if inside or outside a shell
        surrounding the core_mole with radius equal to the cutoff
        threshold
    cut_off: python float number
        The threshold distance in A used to mark atom inside or outside
        the shell
    Return:
    -------
    in_out: python boolean
         True if at least one of check_mol atom distance from core_mole
         is less than the selected cutoff threshold
    """

def sanitizeOEMolecule(molecule):
    """
    This function checks if the molecule has coordinates,
    explicit hydrogens, aromaticity missing and not unique
    atom names. If the molecule does not have coordinates
    a fatal error is raised. If the molecule does not have
    hydrogens or aramatic flags are missing then a copy of
    the molecule is fixed, if missing or not unique atom
    names are found then a copy of the molecule is fixed
    Parameters:
    -----------
    molecule: OEMol
        The molecule to be checked
    Return:
    -------
    mol_copy: OEMol
        A copy of the checked molecule with fixed aromaticity,
        hydrogens and unique atom names if they are missing
    """

def strip_water_ions(in_system):
    """
    This function remove waters and ions molecules
    from the input system
    Parameters:
    ----------
    in_system : oechem.OEMol
        The bio-molecular system to clean
    opt: python dictionary
        The system option
    Output:
    -------
    clean_system : oechem.OEMol
        The cleaned system
    """

def split(complex, ligand_res_name='LIG'):
    """
    This function splits the passed system in protein, ligand,
    water and excipients
    Parameters:
    ----------
    complex : oechem.OEMol
        The bio-molecular complex to split
    ligand_res_name : Python string
        The ligand residue name used to recognize the ligand
    Output:
    -------
    protein : oechem.OEMol
        The split protein
    ligand : oechem.OEMol
        The split ligand
    wat : oechem.OEMol
        The spit water
    other : oechem.OEMol
        The excipients
    """

def select_oemol_atom_idx_by_language(system, mask=''):
    """
    This function selects the atom indexes from the passed oemol molecular complex
    by using  a defined language. The language allows the selection of the ligand,
    protein, waters, ions, cofactors, residue numbers and distance selection. Logic
    operators not, or, and, noh, diff, around can be used to refine the selection
    Parameters
    ----------
    system : OEMol of the bio-molecular complex protein-ligand
        The molecular complex
    mask : python string
        A string used to select atoms. A Backus–Naur Form grammar
        (https://en.wikipedia.org/wiki/Backus–Naur_form) is defined by the python
        module pyparsing.
        The defined grammar tokens are: "ligand", "protein", "ca_protein" ,"water",
        "ions","cofactors" and "resid chain1:res_idx1 chain2:res_idx2 ... res_idxn"
        that respectively define the ligand, the protein, carbon alpha protein atoms,
        water molecules, ions, cofactors and residue numbers. The atom selection can
        be refined by using the following operator tokens:
        "not" = invert selection
        "or" = add selections
        "and" = intersect selections
        "diff" = logic difference between selections
        "noh" = remove hydrogens from the selection
        "around" = select atoms inside the cutoff distance from a given selection
    Returns
    -------
    atom_set : python set
        the select atom indexes
    Notes
    -----
        Example of selection string:
        mask = "ligand or protein"
        mask = "not water or not ions"
        mask = "ligand or protein or cofactors"
        mask = "noh protein"
        mask = "resid A:17 B:12 17 18"
        mask = "protein diff resid A:1"
        mask = "5.0 around protein"
    """
