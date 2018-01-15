from __future__ import division, absolute_imports, print_function
from rdkit import Chem
from simtk.openmm import app, Vec3
from simtk import unit
import itertools
import pyparsing as pyp

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']

def rdkmol_to_openmmTop(mol):
        """
    This function converts an rdkmol to an openmm topology
    The rdkmol coordinates are assumed to be in Angstrom unit
    Parameters:
    -----------
    mol: rdkmol molecule
        The molecule to convert
    Return:
    -------
    topology : OpenMM Topology
        The generated OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the
        generated topology in Angstrom units
    """
    #TODO

def openmmTop_to_rdkmol(topology, positions, verbose = False):
    #TODO

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
