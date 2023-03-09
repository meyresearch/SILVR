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
    recOpts = oedocking.OEMakeReceptorOptions()
    opts = oechem.OESimpleAppOptions(recOpts, "MakeReceptor", oechem.OEFileStringType_DU, oechem.OEFileStringType_DU)
    if oechem.OEConfigureOpts(opts, argv, False) == oechem.OEOptsConfigureStatus_Help:
        return 0
    recOpts.UpdateValues(opts)

    ifs = oechem.oeifstream()
    if not ifs.open(opts.GetInFile()):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % opts.GetInFile())

    ofs = oechem.oeofstream()
    if not ofs.open(opts.GetOutFile()):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % opts.GetOutFile())

    du = oechem.OEDesignUnit()
    while oechem.OEReadDesignUnit(ifs, du):
        if oedocking.OEMakeReceptor(du, recOpts):
            oechem.OEWriteDesignUnit(ofs, du)
        else:
            oechem.OEThrow.Warning("%s: %s" % (du.GetTitle(), "Failed to make receptor"))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
