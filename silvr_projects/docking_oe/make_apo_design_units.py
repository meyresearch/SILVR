#!/usr/bin/env python
# (C) 2018 OpenEye Scientific Software Inc. All rights reserved.
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

#############################################################################
# Script to prepare proteins into design units
#############################################################################
import sys
import os
from openeye import oechem
from openeye import oegrid
from openeye import oespruce


def main(argv=sys.argv):

    if len(argv) < 3 or len(argv) > 5:
        oechem.OEThrow.Usage(
            "%s <infile> <site_residue> [<mtzfile>] [<loopdbfile>]" % argv[0]
        )

    ifs = oechem.oemolistream()
    ifile = argv[1]
    if not ifs.open(ifile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % ifile)

    site_residue = argv[2]

    include_loop = False
    include_ed = False
    if len(argv) > 3:
        if len(argv) == 5 or (len(argv) == 4 and "mtz" in argv[3]):
            edfile = argv[3]
            ed = oegrid.OESkewGrid()
            if not oegrid.OEReadMTZ(edfile, ed, oegrid.OEMTZMapType_Fwt):
                oechem.OEThrow.Fatal(
                    "Unable to read electron density file %s" % edfile
                )  # noqa
            include_ed = True
        if len(argv) == 5:
            loopfile = argv[4]
            include_loop = True
        elif len(argv) == 4 and "mtz" not in argv[3]:
            loopfile = argv[3]
            include_loop = True

    if ifs.GetFormat() not in [oechem.OEFormat_PDB, oechem.OEFormat_CIF]:
        oechem.OEThrow.Fatal("Only works for .pdb or .cif input files")

    ifs.SetFlavor(
        oechem.OEFormat_PDB,
        oechem.OEIFlavor_PDB_Default
        | oechem.OEIFlavor_PDB_DATA
        | oechem.OEIFlavor_PDB_ALTLOC,
    )  # noqa

    mol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ifs, mol):
        oechem.OEThrow.Fatal("Unable to read molecule from %s" % ifile)

    metadata = oespruce.OEStructureMetadata()
    opts = oespruce.OEMakeDesignUnitOptions()
    opts.GetPrepOptions().GetBuildOptions().GetLoopBuilderOptions().SetBuildTails(False)
    if include_loop:
        opts.GetPrepOptions().GetBuildOptions().GetLoopBuilderOptions().SetLoopDBFilename(
            loopfile
        )

    if include_ed:
        design_units = oespruce.OEMakeDesignUnits(mol, ed, metadata, opts, site_residue)
    else:
        design_units = oespruce.OEMakeDesignUnits(mol, metadata, opts, site_residue)

    base_name = os.path.basename(ifile)[:-4] + "_DU_{}.oedu"
    for i, design_unit in enumerate(design_units):
        oechem.OEWriteDesignUnit(base_name.format(i), design_unit)


if __name__ == "__main__":
    sys.exit(main(sys.argv))