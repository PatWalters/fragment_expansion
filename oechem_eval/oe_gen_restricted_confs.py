#!/usr/bin/env python

from docopt import docopt
from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *
from tqdm import tqdm

cmd_str = """Usage:
oe_gen_restricted_confs.py --smi SMILES_FILE --fix FIX_FILE --out OUTPUT_FILE

Options:
--smi SMILES_FILE input SMILES file name
--fix FIX_FILE file_with_fixed piece of the molecule
--out OUTPUT_FILE output file name
"""

cmd_input = docopt(cmd_str)
smiles_file_name = cmd_input.get("--smi")
fix_file_name = cmd_input.get("--fix")
output_file_name = cmd_input.get("--out")

fix_mol = OEGraphMol()
fix_fs = oemolistream(fix_file_name)
OEReadMolecule(fix_fs, fix_mol)

# Relatively quick way to get the number of molecules
moldb = oechem.OEMolDatabase()
moldb.Open(smiles_file_name)
num_mols = moldb.NumMols()

omegaOpts = OEOmegaOptions()
omegaOpts.SetFixMol(fix_mol)
omegaOpts.SetWarts(True)
omega = OEOmega(omegaOpts)

ofs = oemolostream(output_file_name)

smi_fs = oemolistream(smiles_file_name)
num_mols_out = 0
for mol in tqdm(smi_fs.GetOEMols(), total=num_mols):
    if OEGetReasonableProtomer(mol):
        omega.Build(mol)
        if mol.GetDimension() == 3:
            OEWriteMolecule(ofs, mol)
        num_mols_out += 1
print(f"Wrote conformers for {num_mols_out} molecules to {output_file_name}")
