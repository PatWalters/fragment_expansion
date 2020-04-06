from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from tqdm import tqdm
from docopt import docopt
from copy import deepcopy
import Bio.PDB as PDB


# Code borrowed from Joshua Meyers
# https://github.com/JoshuaMeyers/Snippets/blob/master/200405_constrained_conformers.ipynb
# and that code was adapted from Tim Dudgeon
# https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/rdkit/constrained_conf_gen.py
# All I've done is added the commandline wrapper


def duplicate_conformers(m: Chem.rdchem.Mol, new_conf_idx: int, rms_limit: float = 0.5) -> bool:
    rmslist = []
    for i in range(m.GetNumConformers()):
        if i == new_conf_idx:
            continue
        rms = AllChem.GetConformerRMS(m, new_conf_idx, i, prealigned=True)
        rmslist.append(rms)
    return any(i < rms_limit for i in rmslist)


def get_mcs(mol_one: Chem.rdchem.Mol, mol_two: Chem.rdchem.Mol) -> str:
    """Code to find the maximum common substructure between two molecules."""
    return Chem.MolToSmiles(
        Chem.MolFromSmarts(
            rdFMCS.FindMCS([mol_one, mol_two], completeRingsOnly=True, matchValences=True).smartsString
        )
    )


def generate_conformers(mol: Chem.rdchem.Mol,
                        ref_mol: Chem.rdchem.Mol,
                        num_conf: int,
                        ref_smi: str = None,
                        minimum_conf_rms: Optional[float] = None,
                        ) -> List[Chem.rdchem.Mol]:
    # if SMILES to be fixed are not given, assume to the MCS
    if not ref_smi:
        ref_smi = get_mcs(mol, ref_mol)

    # Creating core of reference ligand #
    core_with_wildcards = AllChem.ReplaceSidechains(ref_mol, Chem.MolFromSmiles(ref_smi))
    core1 = AllChem.DeleteSubstructs(core_with_wildcards, Chem.MolFromSmiles('*'))
    core1.UpdatePropertyCache()

    # Add Hs so that conf gen is improved
    mol.RemoveAllConformers()
    outmol = deepcopy(mol)
    mol_wh = Chem.AddHs(mol)

    # Generate conformers with constrained embed
    conf_lst = []
    dup_count = 0
    for i in range(num_conf):
        temp_mol = Chem.Mol(mol_wh)  # copy to avoid inplace changes
        AllChem.ConstrainedEmbed(temp_mol, core1, randomseed=i)
        temp_mol = Chem.RemoveHs(temp_mol)
        conf_idx = outmol.AddConformer(temp_mol.GetConformer(0), assignId=True)
        if minimum_conf_rms is not None:
            if duplicate_conformers(outmol, conf_idx, rms_limit=minimum_conf_rms):
                dup_count += 1
                outmol.RemoveConformer(conf_idx)
    if dup_count:
        pass
    # print(f'removed {dup_count} duplicated conformations')
    return outmol


class ProteinLigandClashFilter:
    def __init__(self, protein_pdbpath: str, distance: float = 1.5):
        parser = PDB.PDBParser(QUIET=True, PERMISSIVE=True)
        s = parser.get_structure('protein', protein_pdbpath)
        self.kd = PDB.NeighborSearch(list(s.get_atoms()))
        self.radius = distance

    def __call__(self, conf: Chem.rdchem.Conformer) -> bool:
        for coord in conf.GetPositions():
            res = self.kd.search(coord, radius=self.radius)
            if len(res):
                return True
        return False


def main():
    cmd_str = """Usage:
    rd_gen_restricted_confs.py --pdb PDB_FILE --smi SMILES_FILE --fix FIX_FILE --out OUTPUT_FILE [--max MAX_CONFS] [--rms RMS]

    Options:
    --pdb PDB_FILE protein pdb file
    --smi SMILES_FILE input SMILES file name
    --fix FIX_FILE file_with_fixed piece of the molecule
    --out OUTPUT_FILE output file name
    --max MAX_CONFS maximum number of conformers to generate - default 25
    --rms RMS RMS cutoff - default 0.01 
    """
    cmd_input = docopt(cmd_str)
    pdb_file_name = cmd_input.get("--pdb")
    smiles_file_name = cmd_input.get("--smi")
    fix_file_name = cmd_input.get("--fix")
    output_file_name = cmd_input.get("--out")
    max_confs = cmd_input.get("--max") or 80
    max_confs = int(max_confs)
    rms = cmd_input.get('--rms') or 0.01
    rms = float(rms)

    ref = Chem.MolFromMolFile(fix_file_name)
    suppl = Chem.SmilesMolSupplier(smiles_file_name, titleLine=False)
    writer = Chem.SDWriter(output_file_name)

    clash_filter = ProteinLigandClashFilter(pdb_file_name, distance=0.5)

    for mol in tqdm(suppl):
        # generate conformers
        out_mol = generate_conformers(mol, ref,
                                      max_confs,
                                      ref_smi=Chem.MolToSmiles(ref),
                                      minimum_conf_rms=rms)

        # remove conformers that clash with the protein
        for conf in out_mol.GetConformers():
            if clash_filter(conf):
                out_mol.RemoveConformer(conf.GetId())

        # write out the surviving conformers
        for conf in out_mol.GetConformers():
            writer.write(out_mol, confId=conf.GetId())


if __name__ == "__main__":
    main()
