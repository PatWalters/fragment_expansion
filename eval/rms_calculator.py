#!/usr/bin/env python

import sys

from openeye import oechem


class RMSCalculator:
    def __init__(self, refmol):
        self.refmol = oechem.OEGraphMol(refmol)
        self.ss = oechem.OESubSearch(oechem.OECreateIsoSmiString(refmol))
        self.ref_match = self.get_match(self.refmol)

    def get_match(self, mol):
        match_list = []
        for mp in self.ss.Match(mol, True):
            for pr in mp.GetAtoms():
                match_list.append(pr.target)
        return match_list

    def calc_rms(self, fitmol):
        fit_match = self.get_match(fitmol)
        match = oechem.OEMatch()
        for rm, fm in zip(self.ref_match, fit_match):
            match.AddPair(rm, fm)
        return oechem.OERMSD(self.refmol, fitmol, match, False)


def main():
    ref_fs = oechem.oemolistream(sys.argv[1])
    fit_fs = oechem.oemolistream(sys.argv[2])

    refmol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ref_fs, refmol)

    rms_calculator = RMSCalculator(refmol)
    for fitmol in fit_fs.GetOEGraphMols():
        print(rms_calculator.calc_rms(fitmol))


if __name__ == "__main__":
    main()
