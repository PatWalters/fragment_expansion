#!/usr/bin/env python

import os

import pandas as pd
from docopt import docopt
from openeye.oechem import *
from openeye.oedocking import *
from openeye.oequacpac import *
from tqdm import tqdm

from bump_check import BumpCheck
from rigid_body_opt import RigidBodyOptimizer
from rms_calculator import RMSCalculator

cmd_str = """Usage:
score_confs.py --prot PROTEIN_FILE --lig LIGAND_FILE --ref REFERENCE_FILE --out OUTPUT_FILE [--bump BUMP_CUTOFF] [--rms RMS_CUTOFF] [--top TOP_OUT]

Options:
--prot PROTEIN_FILE protein file
--lig LIGAND_FILE ligand file to be positioned in the active site
--ref REFERENCE_FILE reference ligand file for comparison
--out OUTPUT_FILE output file name
--bump BUMP_CUTOFF Number of allowed protein-ligand contacts within 2A (default 0)
--rms RMS_CUTOFF only keep structures with with RMS < CUTOFF after minimization (default 1.0 A)
--top TOP_OUT number of structures to write (default 100)
"""


def open_files(cmd_input):
    prot_file_name = cmd_input.get("--prot")
    lig_file_name = cmd_input.get("--lig")
    ref_file_name = cmd_input.get("--ref")
    out_file_name = cmd_input.get("--out")

    prot_fs = oemolistream(prot_file_name)
    lig_fs = oemolistream(lig_file_name)
    ref_fs = oemolistream(ref_file_name)
    out_fs = oemolostream(out_file_name)

    return prot_fs, lig_fs, ref_fs, out_fs


def read_molecules(prot_fs, ref_fs):
    prot_mol = OEGraphMol()
    OEReadMolecule(prot_fs, prot_mol)

    ref_mol = OEGraphMol()
    OEReadMolecule(ref_fs, ref_mol)

    return prot_mol, ref_mol


def build_receptor(prot_fs, prot_mol, ref_mol):
    prot_file_name = prot_fs.GetFileName()
    receptor = OEGraphMol()
    base, _ = os.path.splitext(prot_file_name)
    receptor_file_name = base + "_receptor.oeb"
    if os.path.exists(receptor_file_name):
        print(f"Reading receptor from {receptor_file_name}")
        receptor_fs = oemolistream("receptor.oeb")
        OEReadMolecule(receptor_fs, receptor)
    else:
        print(f"Building receptor")
        OEMakeReceptor(receptor, prot_mol, ref_mol)
        ofs = oemolostream(receptor_file_name)
        OEWriteMolecule(ofs, receptor)
        print(f"Writing receptor to {receptor_file_name}")
    return receptor


def count_molecules(lig_file_name):
    print("Building molecule database")
    moldb = oechem.OEMolDatabase()
    moldb.Open(lig_file_name)
    print("Done")
    num_mols = moldb.NumMols()
    return num_mols


def write_output(score_list, out_fs, num_to_write):
    score_df = pd.DataFrame(score_list, columns=["mol_idx", "name", "score", "rms", "bump", "mol"])
    score_df.sort_values("score", ascending=True, inplace=True)
    score_df.to_csv("scores.csv", index=False)

    num_out = 0
    for i, row in score_df[0:num_to_write].iterrows():
        mol = row.mol
        OESetSDData(mol, "RMS", "%.2f" % row.rms)
        OESetSDData(mol, "Score", "%.2f" % row.score)
        OESetSDData(mol, "Bump", "%d" % row.bump)
        OEWriteMolecule(out_fs, mol)
        num_out += 1
    print(f"Wrote {num_out} structures to {out_fs.GetFileName()}")


def score_molecules(prot_fs, lig_fs, ref_fs, rms_cutoff, bump_cutoff):
    prot_mol, ref_mol = read_molecules(prot_fs, ref_fs)
    receptor = build_receptor(prot_fs, prot_mol, ref_mol)
    score = OEScore(OEScoreType_Chemgauss4)
    bump = BumpCheck(prot_mol)
    rb_opt = RigidBodyOptimizer(prot_mol)
    score.Initialize(receptor)

    score_list = []
    rms_calc = RMSCalculator(ref_mol)
    num_mols = count_molecules(lig_fs.GetFileName())
    for idx, pose in tqdm(enumerate(lig_fs.GetOEGraphMols(), 1), total=num_mols):
        num_bump = bump.count(pose)
        if num_bump <= bump_cutoff:
            opt_result = rb_opt.minimize(pose)
            if opt_result is not None:
                before, after = opt_result
                if after < 0:
                    rms = rms_calc.calc_rms(pose)
                    if rms < rms_cutoff:
                        score_list.append(
                            [idx, pose.GetTitle(), score.ScoreLigand(pose), rms, num_bump, OEGraphMol(pose)])
    return score_list


def main():
    cmd_input = docopt(cmd_str)
    num_to_output = cmd_input.get("--top") or 100
    num_to_output = int(num_to_output)
    rms_cutoff = cmd_input.get("--rms") or 1.0
    rms_cutoff = float(rms_cutoff)
    bump_cutoff = cmd_input.get("--bump") or 0
    bump_cutoff = int(bump_cutoff)
    prot_fs, lig_fs, ref_fs, out_fs = open_files(cmd_input)
    score_list = score_molecules(prot_fs, lig_fs, ref_fs, rms_cutoff, bump_cutoff)
    write_output(score_list, out_fs, num_to_output)


if __name__ == "__main__":
    main()
