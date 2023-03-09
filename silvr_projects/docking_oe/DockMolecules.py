#!/usr/bin/env python
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

import sys

from openeye import oechem
from openeye import oedocking


def main(argv=[__name__]):
    dockOpts = oedocking.OEDockOptions()
    opts = oechem.OERefInputAppOptions(dockOpts, "DockMolecules", oechem.OEFileStringType_Mol3D,
                                       oechem.OEFileStringType_Mol3D, oechem.OEFileStringType_DU, "-receptor")
    if oechem.OEConfigureOpts(opts, argv, False) == oechem.OEOptsConfigureStatus_Help:
        return 0
    dockOpts.UpdateValues(opts)

    ifs = oechem.oemolistream()
    if not ifs.open(opts.GetInFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetInFile())

    rfs = oechem.oeifstream()
    if not rfs.open(opts.GetRefFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetRefFile())

    ofs = oechem.oemolostream()
    if not ofs.open(opts.GetOutFile()):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % opts.GetOutFile())

    du = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(rfs, du):
        oechem.OEThrow.Fatal("Failed to read design unit")
    if not du.HasReceptor():
        oechem.OEThrow.Fatal("Design unit %s does not contain a receptor" % du.GetTitle())

    dock = oedocking.OEDock(dockOpts)
    # @ <SNIPPET-DOCK-MOLECULES-INITIALIZE>
    dock.Initialize(du)
    # @ </SNIPPET-DOCK-MOLECULES-INITIALIZE>

    for mcmol in ifs.GetOEMols():
        print("docking", mcmol.GetTitle())
        dockedMol = oechem.OEGraphMol()
        # @ <SNIPPET-DOCK-MOLECULES-DOCK>
        retCode = dock.DockMultiConformerMolecule(dockedMol, mcmol)
        if (retCode != oedocking.OEDockingReturnCode_Success):
            oechem.OEThrow.Fatal("Docking Failed with error code " + oedocking.OEDockingReturnCodeGetName(retCode))

        # @ </SNIPPET-DOCK-MOLECULES-DOCK>
        sdtag = oedocking.OEDockMethodGetName(dockOpts.GetScoreMethod())
        # @ <SNIPPET-DOCK-MOLECULES-ASSIGN-SCORE>
        oedocking.OESetSDScore(dockedMol, dock, sdtag)
        # @ </SNIPPET-DOCK-MOLECULES-ASSIGN-SCORE>
        # @ <SNIPPET-DOCK-MOLECULES-ANNOTATE>
        dock.AnnotatePose(dockedMol)
        # @ </SNIPPET-DOCK-MOLECULES-ANNOTATE>
        oechem.OEWriteMolecule(ofs, dockedMol)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))