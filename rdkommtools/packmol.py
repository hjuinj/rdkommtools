from __future__ import division, absolute_imports, print_function
"""
openmoltools wrapper for packmol ? https://github.com/choderalab/openmoltools
"""

def oesolvate(solute, density=1.0, padding_distance=10.0,
              distance_between_atoms=2.5,
              solvents='[H]O[H]', molar_fractions='1.0',
              geometry='box', close_solvent=True,
              salt='[Na+], [Cl-]', salt_concentration=0.0,
              neutralize_solute=True, verbose=False, **kargs):
    """
    This function solvates the passed solute in a cubic box or a sphere by using Packmol. Packmol
    creates an initial point for molecular dynamics simulations by packing molecule in defined regions
    of space. For additional info:
    http://www.ime.unicamp.br/~martinez/packmol/home.shtml
    The geometry volume is estimated by the using the padding parameter and the solute size.
    The number of solvent molecules is calculated by using the specified density and volume.
    Solvent molecules are specified as comma separated smiles strings. The molar fractions
    of each solvent molecule are specified in a similar fashion. By default if the solute is
    charged counter ions are added to neutralize it
    Parameters:
    -----------
    solute: OEMol molecule
        The solute to solvate
    density: float
        The solution density in g/ml
    padding_distance: float
        The largest dimension of the solute (along the x, y, or z axis) is determined (in A),
        and a cubic box of size (largest dimension)+2*padding is used
    distance_between_atoms: float
        The minimum distance between atoms in A
    solvents: python string
        A comma separated smiles string of the solvent molecules
    molar_fractions: python string
        A comma separated molar fraction string of the solvent molecules
    close_solvent: boolean
        If True solvent molecules will be placed very close to the solute
    salt: python string
        A comma separated string of the dissociated salt in solution
    salt_concentration: float
        Salt concentration in millimolar
    neutralize_solute: boolean
        If True counter-ions will be added to the solution to neutralize the solute
    Return:
    -------
    oe_mol: OEMol
        The solvated system. If the selected geometry is a box a SD tag with
        name 'box_vector' is attached the output molecule containing
        the system box vectors
    """
