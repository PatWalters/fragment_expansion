#!/usr/bin/env python

import sys

from openeye import oechem


class BumpCheck:
    def __init__(self, prot_mol, cutoff=2.0):
        self.near_nbr = oechem.OENearestNbrs(prot_mol, cutoff)
        self.cutoff = cutoff

    def count(self, lig_mol):
        bump_count = 0
        for nb in self.near_nbr.GetNbrs(lig_mol):
            if (not nb.GetBgn().IsHydrogen()) and (not nb.GetEnd().IsHydrogen()) and nb.GetDist() <= self.cutoff:
                bump_count += 1
        return bump_count


if __name__ == "__main__":
    PROT = sys.argv[1]
    LIG = sys.argv[2]
    prot = oechem.OEMol()
    oechem.OEReadPDBFile(oechem.oemolistream(PROT), prot)

    bump_check = BumpCheck(prot)
    for lig_mol in oechem.oemolistream(LIG).GetOEGraphMols():
        print(bump_check.count(lig_mol))
