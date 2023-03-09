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

#############################################################################
# prepare a molecule: alts, hydrogens, split ligand
#############################################################################
import sys
from openeye import oechem


def WaterProcess(processName):
    if processName == "fullsearch":
        return oechem.OEPlaceHydrogensWaterProcessing_FullSearch
    elif processName == "focused":
        return oechem.OEPlaceHydrogensWaterProcessing_Focused
    return oechem.OEPlaceHydrogensWaterProcessing_Ignore


def main(argv=[__name__]):

    itf = oechem.OEInterface(InterfaceData)
    oechem.OEConfigureSplitMolComplexOptions(itf,
                                             oechem.OESplitMolComplexSetup_All &
                                             ~ (oechem.OESplitMolComplexSetup_CovBondTreatment |
                                                oechem.OESplitMolComplexSetup_CovCofactor))

    if not oechem.OEParseCommandLine(itf, argv):
        oechem.OEThrow.Fatal("Unable to interpret command line!")

    verbose = itf.GetBool("-verbose")
    if verbose:
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Verbose)

    altProcess = itf.GetString("-alts")
    keepAlts = (altProcess != "a")
    highestOcc = (altProcess == "occupancy")
    compareAlts = (altProcess == "compare")

    siteNum = itf.GetUnsignedInt("-bindingsitenum")
    allSites = (siteNum == 0)
    otherModel = (itf.GetUnsignedInt("-modelnum") != 1)

    placeHyd = itf.GetBool("-placehydrogens")
    splitlig = itf.HasString("-ligout")

    watProcessName = itf.GetString("-waterprocessing")
    waterProcess = WaterProcess(watProcessName)

    standardize = itf.GetBool("-standardizehyd")
    badclash = itf.GetDouble("-clashcutoff")
    flipbias = itf.GetDouble("-flipbias")
    maxStates = itf.GetDouble("-maxsubstates")

    flavor = oechem.OEIFlavor_PDB_Default | oechem.OEIFlavor_PDB_DATA
    if keepAlts:
        flavor = flavor | oechem.OEIFlavor_PDB_ALTLOC
    if otherModel:
        flavor = flavor & ~ oechem.OEIFlavor_PDB_ENDM

    ims = oechem.oemolistream()
    ims.SetFlavor(oechem.OEFormat_PDB, flavor)

    inputFile = itf.GetString("-in")
    if not ims.open(inputFile):
        oechem.OEThrow.Fatal("Unable to open %s for reading." % inputFile)

    if not oechem.OEIs3DFormat(ims.GetFormat()):
        oechem.OEThrow.Fatal("%s is not in a 3D format." % inputFile)

    inftype = oechem.OEGetFileType(oechem.OEGetFileExtension(inputFile))
    if (inftype == oechem.OEFormat_PDB) and not keepAlts:
        oechem.OEThrow.Verbose("Default processing of alt locations (keep just 'A' and ' ').")

    sopt = oechem.OESplitMolComplexOptions()
    oechem.OESetupSplitMolComplexOptions(sopt, itf)

    inmol = oechem.OEGraphMol()
    if not oechem.OEReadMolecule(ims, inmol):
        oechem.OEThrow.Fatal("Unable to read %s." % inputFile)

    ims.close()

    if (inmol.NumAtoms() == 0):
        oechem.OEThrow.Fatal("Input molecule %s contains no atoms." % inputFile)

    if inmol.GetTitle() == "":
        inmol.SetTitle("input mol")

    oechem.OEThrow.Verbose("Processing %s." % inmol.GetTitle())

    if not oechem.OEHasResidues(inmol):
        oechem.OEPerceiveResidues(inmol, oechem.OEPreserveResInfo_All)

    if highestOcc or compareAlts:
        alf = oechem.OEAltLocationFactory(inmol)
        if alf.GetGroupCount() != 0:
            if highestOcc:
                oechem.OEThrow.Verbose("Dropping alternate locations from input.")
                alf.MakePrimaryAltMol(inmol)
            elif compareAlts:
                oechem.OEThrow.Verbose("Fixing alternate location issues.")
                inmol = alf.GetSourceMol()

    outmol = oechem.OEGraphMol()
    if allSites:
        outmol = inmol
    else:
        oechem.OEThrow.Verbose("Splitting out selected complex.")

        soptSiteSel = oechem.OESplitMolComplexOptions(sopt)
        soptSiteSel.SetSplitCovalent(False)  # do any cov lig splitting later

        frags = oechem.OEAtomBondSetVector()
        if not oechem.OEGetMolComplexFragments(frags, inmol, soptSiteSel):
            oechem.OEThrow.Fatal("Unable to fragment %s." % inmol.GetTitle())

        howManySites = oechem.OECountMolComplexSites(frags)
        if howManySites < siteNum:
            oechem.OEThrow.Warning(("Binding site count (%d) " +
                                   "less than requested site (%d) in %s.") %
                                   (howManySites, siteNum, inmol.GetTitle()))
            exit(0)

        if not oechem.OECombineMolComplexFragments(outmol, frags, soptSiteSel):
            oechem.OEThrow.Fatal("Unable to collect fragments from %s." % inmol.GetTitle())

        if (outmol.NumAtoms() == 0):
            oechem.OEThrow.Fatal("No fragments selected from %s." % inmol.GetTitle())

    if placeHyd:
        oechem.OEThrow.Verbose("Adding hydrogens to complex.")

        hopt = oechem.OEPlaceHydrogensOptions()
        hopt.SetAltsMustBeCompatible(compareAlts)
        hopt.SetStandardizeBondLen(standardize)
        hopt.SetWaterProcessing(waterProcess)
        hopt.SetBadClashOverlapDistance(badclash)
        hopt.SetFlipBiasScale(flipbias)
        hopt.SetMaxSubstateCutoff(maxStates)

        if verbose:
            details = oechem.OEPlaceHydrogensDetails()
            if not oechem.OEPlaceHydrogens(outmol, details, hopt):
                oechem.OEThrow.Fatal("Unable to place hydrogens and get details on %s."
                                     % inmol.GetTitle())
            oechem.OEThrow.Verbose(details.Describe())
        else:
            if not oechem.OEPlaceHydrogens(outmol, hopt):
                oechem.OEThrow.Fatal("Unable to place hydrogens on %s." % inmol.GetTitle())

    oms1 = oechem.oemolostream()
    cplxFile = itf.GetString("-cplxout")
    if not oms1.open(cplxFile):
        oechem.OEThrow.Fatal("Unable to open %s for writing." % cplxFile)

    if splitlig:
        oechem.OEThrow.Verbose("Splitting ligand from complex.")

        returnAllSites = 0
        oechem.OESetupSplitMolComplexOptions(sopt, itf, returnAllSites)

        frags = oechem.OEAtomBondSetVector()
        if not oechem.OEGetMolComplexFragments(frags, outmol, sopt):
            oechem.OEThrow.Fatal("Unable to fragment complex from %s." % inmol.GetTitle())

        lfilter = sopt.GetLigandFilter()

        protComplex = oechem.OEGraphMol()
        if not oechem.OECombineMolComplexFragments(protComplex,
                                                   frags,
                                                   sopt,
                                                   oechem.OENotRoleSet(lfilter)):
            oechem.OEThrow.Fatal("Unable to collect complex from %s." % inmol.GetTitle())

        if (protComplex.NumAtoms() == 0):
            oechem.OEThrow.Warning("No complex identified in %s." % inmol.GetTitle())
        else:
            oechem.OEWriteMolecule(oms1, protComplex)

        lig = oechem.OEGraphMol()
        if not oechem.OECombineMolComplexFragments(lig, frags, sopt, lfilter):
            oechem.OEThrow.Fatal("Unable to collect ligand from %s." % inmol.GetTitle())

        if (lig.NumAtoms() == 0):
            oechem.OEThrow.Warning("No ligand identified in %s." % inmol.GetTitle())
        else:
            oms2 = oechem.oemolostream()
            if splitlig:
                ligFile = itf.GetString("-ligout")
                if not oms2.open(ligFile):
                    oechem.OEThrow.Fatal("Unable to open %s for writing." % ligFile)

            oechem.OEThrow.Verbose("Ligand: %s" % lig.GetTitle())
            oechem.OEWriteMolecule(oms2, lig)
            oms2.close()
    else:
        oechem.OEWriteMolecule(oms1, outmol)

    oms1.close()

#############################################################################
# INTERFACE
#############################################################################


InterfaceData = '''
!BRIEF proteinprep.py [-options] <inmol> [<outcplx> [<outlig>]]

!CATEGORY "input/output options :" 1
   !PARAMETER -in 1
      !ALIAS -i
      !TYPE string
      !BRIEF Input molecule filename (must have 3D coordinates)
      !SIMPLE true
      !REQUIRED true
      !KEYLESS 1
   !END

   !PARAMETER -cplxout 2
      !ALIAS -p
      !TYPE string
      !DEFAULT proteinprep.oeb.gz
      !BRIEF Output complex filename
      !SIMPLE true
      !REQUIRED false
      !KEYLESS 2
   !END

   !PARAMETER -ligout 3
      !ALIAS -l
      !TYPE string
      !BRIEF Output ligand filename
      !SIMPLE true
      !REQUIRED false
      !KEYLESS 3
   !END
!END

!CATEGORY "Calculation options :" 2
    !PARAMETER -alts 1
       !TYPE string
       !LEGAL_VALUE occupancy
       !LEGAL_VALUE a
       !LEGAL_VALUE ignore
       !LEGAL_VALUE compare
       !DEFAULT occupancy
       !BRIEF Alternate location atom handling (affects atom:atom interactions)
       !SIMPLE true
       !REQUIRED false
       !DETAIL
         occupancy - keep just the highest average occupancy for each alt group
         a - keep only loc code A (and blank)
         ignore - assume alts already selected appropriately
         compare - keep all alts but only interact if same loc code (or blank)
    !END

    !PARAMETER -placehydrogens 2
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, hydrogens will not be added
       !SIMPLE true
       !REQUIRED false
    !END

    !PARAMETER -waterprocessing 3
       !TYPE string
       !LEGAL_VALUE ignore
       !LEGAL_VALUE focused
       !LEGAL_VALUE fullsearch
       !DEFAULT fullsearch
       !BRIEF How waters are processed
       !SIMPLE true
       !REQUIRED false
       !DETAIL
         ignore - leave water hydrogens in a random orientation
         focused - search orientations based on neighboring polar groups
         fullsearch - do an extensive search of water orientations
    !END

    !PARAMETER -standardizehyd 4
       !ALIAS -stdhyd
       !TYPE bool
       !DEFAULT true
       !BRIEF If false, bonds for hydrogens are not adjusted to standard lengths
       !SIMPLE false
       !REQUIRED false
    !END

    !PARAMETER -clashcutoff 5
        !TYPE double
        !DEFAULT 0.4
        !BRIEF Van der Waals overlap (in Angstroms) defined to be a bad clash
        !SIMPLE false
        !REQUIRED false
    !END

    !PARAMETER -flipbias 6
        !TYPE double
        !DEFAULT 1.0
        !BRIEF Scale factor for the bias against flipping sidechains such as HIS
        !SIMPLE false
        !REQUIRED false
    !END

    !PARAMETER -maxsubstates 7
        !TYPE double
        !DEFAULT 1.0e8
        !BRIEF Maximum number of substates in a single step of hydrogen placement optimization
        !SIMPLE false
        !REQUIRED false
    !END
!END

!CATEGORY "Display options :" 9
   !PARAMETER -verbose 1
      !ALIAS -v
      !TYPE bool
      !DEFAULT false
      !BRIEF Display more information about the process
      !SIMPLE true
      !REQUIRED false
   !END
!END
'''

if __name__ == "__main__":
    sys.exit(main(sys.argv))
