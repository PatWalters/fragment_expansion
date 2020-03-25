#!/usr/bin/env python

import sys
import pandas as pd
from rdkit import Chem
from tqdm.notebook import tqdm

class REOS:
    def __init__(self,rules_file):
        self.df = pd.read_csv(rules_file)
        self.df['pat'] = [Chem.MolFromSmarts(x) for x in tqdm(self.df.smarts,desc="Reading rules")]

    def eval_mol(self,mol):
        for pat,max_val in self.df[['pat','max']].values:
            if len(mol.GetSubstructMatches(pat)) > max_val:
                return False
        return True

    def eval_smiles(self,smiles):
        mol = Chem.MolFromSmiles(smiles)
        return self.eval_mol(mol)
        

if __name__ == "__main__":
    reos = REOS("pw_rules.csv")
    suppl = Chem.SmilesMolSupplier(sys.argv[1],titleLine=False)
    for mol in suppl:
        print(reos.eval_mol(mol))
