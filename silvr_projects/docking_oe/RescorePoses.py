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
    itf = oechem.OEInterface(InterfaceData)
    # @ <SNIPPET-RESCORE-POSES-CONFIGURE>
    oedocking.OEScoreTypeConfigure(itf, "-score")
    # @ </SNIPPET-RESCORE-POSES-CONFIGURE>
    if not oechem.OEParseCommandLine(itf, argv):
        return 1

    receptor = oechem.OEDesignUnit()
    if not oechem.OEReadDesignUnit(itf.GetString("-receptor"), receptor):
        oechem.OEThrow.Fatal("Unable to read receptor")
    imstr = oechem.oemolistream()
    if not imstr.open(itf.GetString("-in")):
        oechem.OEThrow.Fatal("Unable to open input file of ligands")
    omstr = oechem.oemolostream()
    if not omstr.open(itf.GetString("-out")):
        oechem.OEThrow.Fatal("Unable to open out file for rescored ligands")
        
        
    scoreType = oedocking.OEScoreTypeGetValue(itf, "-score")
    print("scoreType: ", scoreType)

    score = oedocking.OEScore(scoreType)
    print("score: ", score)

    score.Initialize(receptor)
    

    for ligand in imstr.GetOEMols():
        if itf.GetBool("-optimize"):
            # @ <SNIPPET-RESCORE-POSES-OPTIMIZE>
            score.SystematicSolidBodyOptimize(ligand)
            # @ </SNIPPET-RESCORE-POSES-OPTIMIZE>
        # @ <SNIPPET-RESCORE-POSES-ANNOTATE>
        score.AnnotatePose(ligand)
        # @ </SNIPPET-RESCORE-POSES-ANNOTATE>
        sdtag = score.GetName()
        # @ <SNIPPET-RESCORE-POSES-ASSIGN-SCORE>
        oedocking.OESetSDScore(ligand, score, sdtag)
        # @ </SNIPPET-RESCORE-POSES-ASSIGN-SCORE>
        # @ <SNIPPET-RESCORE-POSES-SCORE-SORTING>
        oechem.OESortConfsBySDTag(ligand, sdtag, score.GetHighScoresAreBetter())
        # @ </SNIPPET-RESCORE-POSES-SCORE-SORTING>
        oechem.OEWriteMolecule(omstr, ligand)

    return 0


InterfaceData = """
!PARAMETER -receptor
  !ALIAS -rec
  !TYPE string
  !REQUIRED true
  !LEGAL_VALUE *.oedu
  !BRIEF A receptor file the poses pass to the -in flag will be scored against
!END

!PARAMETER -in
  !TYPE string
  !REQUIRED true
  !BRIEF Input molecule file with poses to rescore
!END

!PARAMETER -out
  !TYPE string
  !REQUIRED true
  !BRIEF Rescored molecules will be written to this file
!END

!PARAMETER -optimize
  !ALIAS -opt
  !TYPE bool
  !DEFAULT false
  !BRIEF Optimize molecules before rescoring
!END
"""


if __name__ == "__main__":
    sys.exit(main(sys.argv))
