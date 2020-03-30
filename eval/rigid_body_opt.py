#!/usr/bin/env python

import sys

from openeye import oechem
from openeye import oeff


class RigidBodyOptimizer:
    def __init__(self,prot_mol):
        oechem.OEAddExplicitHydrogens(prot_mol)
        self.mmff = oeff.OEMMFFAmber(prot_mol)
        self.mmff.PrepMol(prot_mol)

    def minimize(self,mol):
        oechem.OEAddExplicitHydrogens(mol)
        if (not self.mmff.PrepMol(mol)) or (not self.mmff.Setup(mol)):
            oechem.OEThrow.Warning("Unable to process molecule: title = '%s'" % mol.GetTitle())
            return None

        adaptor = oeff.OEQuatAdaptor(self.mmff, False, False)
        if not adaptor.Setup(mol):
            oechem.OEThrow.Warning("Unable to process subset for molecule: title = '%s'"
                                   % mol.GetTitle())
            return None

        vecCoords = oechem.OEDoubleArray(3*mol.GetMaxAtomIdx())
        mol.GetCoords(vecCoords)

        vecX = oechem.OEDoubleArray(adaptor.NumVar())
        adaptor.GetVar(vecX, vecCoords)

        initial_energy = adaptor(vecX)

        optimizer = oeff.OEBFGSOpt()
        final_energy = optimizer(adaptor, vecX, vecX)

        adaptor.AdaptVar(vecCoords, vecX)
        mol.SetCoords(vecCoords)

        return initial_energy, final_energy


def main():
    prot_fs = oechem.oemolistream(sys.argv[1])
    lig_fs =  oechem.oemolistream(sys.argv[2])
    ofs = oechem.oemolostream(sys.argv[3])

    prot_mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(prot_fs,prot_mol)

    rigid_body_opt = RigidBodyOptimizer(prot_mol)

    for lig_mol in lig_fs.GetOEMols():
        res = rigid_body_opt.minimize(lig_mol)
        print(lig_mol.GetTitle(),res)
        oechem.OEWriteMolecule(ofs,lig_mol)


if __name__ == "__main__":
    main()

