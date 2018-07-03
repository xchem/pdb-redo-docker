#!/bin/tcsh -f

# pdb_redo.csh: a general method for optimising macromolecular crystallographic structures.
#
# PDB_REDO installation directory. EDIT THIS!
setenv BASE /usr/local/pdb-redo
setenv PATH /ccp4/include:/ccp4/bin:/ccp4/lib:/ccp4/lib/data:/ccp4/etc:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
# This script was created by Robbie P. Joosten, Bart van Beusekom, and Wouter Touw
# Correspondence to r.joosten@nki.nl or robbie_joosten@hotmail.com
#
# Reference: If you publish results (directly or indirectly) obtained by using this protocol, please refer to (any of)
# these publications:
# 1) Robbie P. Joosten, Gert Vriend: "PDB improvement starts with data deposition" Science, 317, p. 195-196 (2007)
# 2) Robbie P. Joosten, Thomas Womack, Gert Vriend and Gerard Bricogne: "Re-refinement fromdeposited X-ray data can
#    deliver improved models for most PDB entries"  Acta Cryst. D65, p. 176-185 (2009)
# 3) Robbie P. Joosten, Jean Salzemann, Vincent Bloch, Heinz Stockinger, Ann-Charlott Berglund, Christophe Blanchet, Erik
#    Bongcam-Rudloff, Christophe Combet, Ana L. Da Costa, Gilbert Deleage, Matteo Diarena, Roberto Fabbretti, Geraldine
#    Fettahi, Volker Flegel, Andreas Gisel, Vinod Kasam, Timo Kervinen, Eija Korpelainen, Kimmo Mattila, Marco Pagni,
#    Matthieu Reichstadt, Vincent Breton, Ian J. Tickle, Gert Vriend: "PDB_REDO: automated re-refinement of X-ray
#    structure models in the PDB" J. Appl. Cryst., 42, p. 376-384 (2009)
# 4) Robbie P. Joosten, Tim A.H. te Beek, Elmar Krieger, Maarten Hekkelman, Rob W.W. Hooft, Reinhard Schneider, Chris
#    Sander, Gert Vriend: "A series of PDB related databases for everyday needs" Nucl. Acids Res., 39, p. D411-D419 (2011)
# 5) Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, Anastassis Perrakis: "Automatic rebuilding and
#    optimisation of crystallographic structures in the Protein Data Bank" Bioinformatics, 27, p. 3392-3398 (2011)
# 6) Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastass Perrakis: "PDB_REDO: constructive validation, more
#    than just looking for errors" Acta Cryst. D68, p. 484-496 (2012)
#

# Usage notice:

if ($#argv == 0 ) then
  echo "How to use PDB_REDO:"
  echo " "
  echo "Optimising your structure model (locally):"
  echo "'$0 --local --pdbin=a_pdb_file --hklin=an_mmCIF_file (or --mtzin=an_mtz_file) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--seqin=fasta_sequence_file) --dirout=place_for_output (--flags)'"
  echo " "
  echo "Optimising your structure model (as a PDB_REDO webserver):"
  echo "'$0 --server --pdbin=a_pdb_file --hklin=an_mmCIF_file (or --mtzin=an_mtz_file) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--seqin=fasta_sequence_file) --dirout=place_for_output (--flags)'"
  echo " "
  echo "Optimising an entry from the Protein Data Bank (PDB):"
  echo "'$0 [PDBid] (--download) (--tlsin=tls_file) (--homin=pdb_file) (--restin=cif_restraint_file) (--extin=external_restraints) (--dirout=place_for_output) (--flags)'"
  echo " "
  echo "Please note:"
  echo "-The --restin option takes a single restraint file for Refmac in mmCIF format."
  echo " It overrides all the automatic restraint generation."
  echo "-The --extin option takes a single external restraint file in Refmac format."
  echo " The external restraints are added to those generated within PDB_REDO."
  echo "-The --seqin option takes a sequence file in fasta format."
  echo " The sequence is cross checked with the (ATOM and SEQRES) sequence in the pdb file."
  echo "-The --tlsin option takes an extra TLS group definition in Refmac format."
  echo " You can give the --tlsin option multiple times to test more TLS group definitions."
  echo "-The --tlsin option takes an extra TLS group definition in Refmac format."
  echo " You can give the --tlsin option multiple times to test more TLS group definitions."  
  echo "-You can add a comment to the output HTML file: --comment='This is a comment'"
  echo "-The --download flag lets the required data from be downloaded from PDBe. Without"
  echo " this flag the required data is taken from a local copy of the PDB."
  echo " "
  echo "The debug flags do the following:"
  echo "--nproc=VAL   : the maximum number of CPU cores you want to use (the default is 1);"
  echo "                values greater than 7 or the actual number of CPU cores on your"
  echo "                system are not recommended"
  echo "--nohyd       : do not add hydrogens (in riding postions) during refinement"
  echo "--legacy      : for legacy PDB entries. R-factor is not checked and the number of "
  echo "                refinement cycles is increased (a lot)"
  echo "--notls       : no TLS refinement is performed"
  echo "--notlsupdate : use TLS, but do not update the tensors in the final refinement"
  echo "--noncs       : no NCS restraints are applied"
  echo "--nojelly     : switch off jelly body refinement"
  echo "--notwin      : no detwinning is performed"
  echo "--newmodel    : always take an updated model from the re-refinement for the rebuilding"
  echo "                steps. Only use this option when all else fails"
  echo "--tighter     : try tighter restraints than usual (use '--tighter --tighter' for even"
  echo "                tighter restraints)"
  echo "--looser      : try looser restraints than usual (use '--looser --looser' for even"
  echo "                looser restraints)"
  echo "--nopepflip   : no peptide flips are performed"
  echo "--noscbuild   : side chains will not be rebuilt"
  echo "--nocentrifuge: waters with poor density will not be deleted"
  echo "--norebuild   : all rebuilding steps are skipped"
  echo "--noanomalous : ignore all anomalous data if Fmean or Imean are available"
  echo "--lowmem      : reduce memory consumption for massive models that otherwise stop Refmac"
  echo "--maxres=VAL  : cut the resolution to VAL Angstrom"
  echo "--fewrefs     : deals with very small data sets by switching off R-free set sanity checks"
  echo "--crossval    : performs (very lengthy) k-fold cross validation on the final results"
  echo "--intens      : force the use of the intensities from the reflection file"
  echo "--noocc       : do not refine occupancies"
  echo "--notruncate  : do not use truncate to convert intensities to amplitudes"
  echo "--nosigma     : do not use sigF or sigI for scaling etc."
  echo "--nometalrest : do not generate special metal restraints"
  echo "--nohomology  : do not use homology-based restraints"
  echo "--homology    : force homology-based restraints"
  echo "--hbondrest   : use hydrogen bond restraints"
  echo "--libg        : use DNA/RNA restraints from LibG"
  echo " "
  echo "Debug flags only for working on local PDB entries:"
  echo "--nopdb       : any existing PDB file in the working directory is not replaced"
  echo "--nosf        : any existing reflection file in the working directory is not replaced"
  echo "--mtzdb       : the structure factor files are stored in MTZ format"
  echo " "
  echo "Citing PDB_REDO:"
  echo "-Robbie P. Joosten, Krista Joosten, Serge X. Cohen, Gert Vriend, Anastassis Perrakis: "
  echo " Automatic rebuilding and optimization of crystallographic structures in the Protein Data Bank"
  echo " Bioinformatics, 27, p. 3392-3398 (2011)"
  echo "-Robbie P. Joosten, Krista Joosten, Garib N. Murshudov, Anastassis Perrakis:"
  echo " PDB_REDO: constructive validation, more than just looking for errors"
  echo " Acta Cryst. D68, p. 484-496 (2012)"
  exit(0)

endif
#
# Required software:
# - PDB_REDO       Download it from http://www.cmbi.ru.nl/pdb_redo
# - CCP4 package   Download it from http://www.ccp4.ac.uk
# - WHAT_CHECK     Download it from http://swift.cmbi.ru.nl/gv/whatcheck
# - pdb-care       See http://www.glycosciences.de for more information
#
# Optional software
# - BlastP         Download it from or apt install ncbi-blast+
#                  Needed for homology-based functionality
# - YASARA         Download it from http://www.yasara.org
#                  You need at least YASARA model for pictures describing atomic shifts and TLS grouping and
#                  YASARA dynamics for ligand validation
# - FoldX          Download it from http://foldxsuite.crg.es
#                  Needed for protein stability analysis
# - R (statistics) Download it from http://www.r-project.org or apt install r-base
#                  Needed for density fit analysis and for boxplots on relative model quality
#
####################################################### Change log #######################################################
set VERSION = '7.08' #PDB_REDO version

# Version 7.08:
# - Added restraint geration with KRAB. Use the --krab switch to activate it.
# - Changed the setting for restraint generation in Refmac to 'connectivity YES'.
# - Simplified the way the version-tracking file is made.
#
# Version 7.07:
# - If the R-factor cannot be reproduced in databank mode, but a previous PDB-REDO entry exists, the 0-cycle coordinates 
#   are downloaded and used to try to reproduce the R-factors.
# - The hidden --no rb switch switches off rigid-body refinement, even in legacy mode.
#
# Version 7.06:
# - More elegant overwriting of databank entries. This does mean that old data in the output directory is removed.
# - The PDBe.json files are now always created.
# - The index.html files are no longer created.
# - WHYNOT and DEBUG messages for server runs go to a seperate WHYNOT file.  
#
# Version 7.05:
# - Updated data for percentile calculation and boxplots.
# - Now also calculating percentiles for coarse packing.
# - CB atoms of GLY are removed automatically.
#
# Version 7.04:
# - REMARK 350 records (biological assembly) are now retained if available.
# - Workaround for cases with alternate residues that get mangled by the rebuilding tools.
# - Workaround for cases with a rediculous amount of TLS groups.
# - Better handling of carbohydrate LINKs with alternates.
#
# Version 7.03:
# - With '--homin=homologous_pdb_file' users can add additional homologous structures that are not yet in the PDB to 
#   improve the homology restraints. 
# - This option automatically triggers homology restraints.
# - Bugfix in writing out the refmac script for low resolution cases.
# - The DNA/RNA restraints from LibG are no longer a hidden option.
# - Fix importing external restraints.
#
# Version 7.02:
# - Added symmetry support for Zen. This required generation of a temporary mtz file early in PDB-REDO.
# - More elegant handling of TLS groups not working in TLSanl. 
# - Experimental phases without HL coefficients or FOMs are supplemented with default FOMs.
# - If the NCS alignment failes due to terminal insertion codes. The residues are now automatically renumbered.
# - Bugfix for cases that have problems with ctruncate and have relatively small test sets.
#
# Version 7.01:
# - Changed the slider values for PDBe to make sure that the step from neutral to positive happens for values greater than
#   2.58.
# - More elegant handling of unstable TLS refinement.
# - More elegant handling of cases where no TLS file could be made.
# - Better handling of cases where Refmac refuses to make anomalous maps.
# - More elegant handling of cases where LibG does not make all types of restraints. 
# - Bug fix in the way Phaser errors are intercepted.
# - Fix for reading logfiles of backbone-only models. 
#
# Version 7.00:
# - rmsZ values for (homology-based) H-bond restraints are now reported. 
# - Better handling for missing high resolution data.
# - Changed the way data.txt is written to avoid race conditions.
#
# Version 6.29:
# - The output PDB files now also have SEQRES records.
# - Added flipper to do DEFY standardisation flips. This replaces the WHAT_CHECK run on the not-quite final pdb file.
# - Flipper is also run before extractor to avoid potential harmful DEFY flips by other programs.
# - Workaround for ctruncate problems if there are anomalous intensities.
#
# Version 6.28:
# - Input MTZ files with many systematic absent reflection now no longer cause infinite loops.
# - Fix for twinning detection by PHASER.
# - Specified NOHARVEST for the REFMAC jobs.
# - Fix for new version of WHAT_CHECK.
#
# Version 6.27:
# - Homology restraints are now automatically used in the third lowest resolution category.
# - Added the --nohomology flag to not use homology restraints.
# - Added H-bond satisfaction percentile for the server.
# - Writing out percentiles to data.txt.
# - Removed the link to EDS.
# - More robust treatment of alternative directory structures and compression of local PDB/reflection files.
# - The besttls files are now gzipped in the databank.
# - Missing percentiles are no longer in <em> tags on the server.
#
# Version 6.26:
# - Homology restraints are now also used in the second lowest resolution category.
# - Fix to deal with = characters in the file names.
# - Carbohydrate residues are now only renamed if they are part of N-glycans.
# - PDB-care is run differently based on the presence or absence of CONECT records.
#
# Version 6.25:
# - Added a twin test from PHASER to reduce twinning false positives.
# - Fixed bug in handeling corrupt restraint files.
# - SEGID records (not part of the official PDB format) are now deleted.
# - Added the --mtzdb flag.
#
# Version 6.24:
# - The configuration of PDB_REDO is moved to a separate file which makes updating easier.
# - Changed the geometric weight for mini-rsr to 10.00 based on feedback from Paul Emsley.
# - Users can give a sequence file.
# - Now using Refmac's default treatment of sugars.
# - Stability fix for parsing WHAT_CHECK output.
# - Better handling of incorrect external restraint files.
#
# Version 6.23:
# - Fix in the map handeling for pepflip.
# - Improved error message for atom naming problems.
# - Bug fix for PHASER output detection.
# - Kollumer now ignores anomalous data with very low completeness.
# - Homology restraints are now automatically used in the lowest resolution category.
# - If the resolution cut-off is changed the solvent mask parameters are re-optimised.
# - Bugfix in YASARA atom shift scene generation.
#
# Version 6.22:
# - A JSON version of data.txt is now also written.
# - Started using WHAT_CHECK 14.
# - Stability fix for homology restraints.
# - Better handling of rediculously high (total) B-factors after TLS refinement.
#
# Version 6.21:
# - Updated the information in the data.txt that allows regeneration of all relevant warnings in the index.html files.
# - The entry creation date is now also stored in data.txt.
# - Stabitility fix in the real-space map validation.
# - Fixed support for users without FoldX.
# - Substantial reduction in messages to debug.txt.
#
# Version 6.20:
# - Added workaround for Refmac deleting LINK records after rigid-body refinement.
# - If ctruncate fails or writes out an unusable reflection file the I to F conversion is performed by cif2cif.
# - Models with strict NCS are now also automatically run in databank mode.
#
# Version 6.19:
# - Ligand validation in YASARA is now also parallelised.
# - Fixed bug in the residue name extraction for validated ligands.
# - The wavelength is now read form the PDB header if no value is given in the reflection file.
# - Stability fix for k-fold cross validation.
# - Switched to the CCP4 version of syminfo.lib for the rebuilding tools.
#
# Version 6.18:
# - Switched back to the mini-rsr (now named coot-mini-rsr) distributed with CCP4.
# - At low resolution the B-factor model selection is slightly less conservative.
# - A debug statement is written for legacy entries with calculated R-free > 0.50.
# - PDB_REDO is stopped if the resolution of the reflection data is higher than 0.30A, which is an indication reflection
#   data problems.
# - Stability fix in the real-space validation.
#
# Version 6.17:
# - Anomalous data is now used to make anomalous maps, but not for refinement.
# - Generating anomalous maps is not compatible with twinning or phased refinement.
# - Anomalous data can be ignored with the '--noanomalous' flag.
# - Small improvements to how the MTZ file is made.
#
# Version 6.16:
# - The completeness of the data is now reported in data.txt.
# - At very low completeness (less than 50%) and omit map is created with COMIT. This map is used for water deletion,
#   side-chain rebuilding and selection of peptide flipping candidates. The map is of too low quality for mini-rsr.
# - Made centrifuge a bit less aggressive in water deletion.
# - The new program rotacompare replaces YASARA for rotamer and peptide analysis.
# - External restraint files are now checked to see whether they are not regular mmCIF restraint files.
#
# Version 6.15:
# - Change in pepflip that should improve the initial filtering of flipping candidates.
# - Changed resolution cut-off for peptide flipping to 3.3A (from 3.5A).
# - Switched to a new rotamer dictionary for SideAide, based on the Top8000/Top7200.
# - Reporting of missing wavelengths is now limited to entries from 2007 onwards.
#
# Version 6.14:
# - A JSON-formatted datafile for PDBe is now written in databank mode.
# - A debug message is written when the wavelength is not reported in the reflection file.
# - Residues named WAT are now excluded from ligand validation, UNL is now included.
# - The output from the rebuilding tools is now slightly less verbose.
# - Stability fix in making the dRSCC plots.
#
# Version 6.13:
# - The solvent percentage is now also given for structures with strict NCS.
# - Hydrogen bond restraints now also apply to side chains.
# - Added optional homology based restraints, this requires a local installation of BlastP.
# - Homology-based restraints currently only work in databank mode.
# - Activate homology-based restraints with --homology.
#
# Version 6.12:
# - LINKs are now added based on the output of PDB-care.
#
# Version 6.11:
# - CIF2CIF now keeps phase information from the input reflection file.
# - Several related stability fixes to cif2cif.
# - Started using CCP4 7.0.
#
# Version 6.10:
# - First attempt to deal with electron diffraction data. It is not very sophisticated yet.
# - The type of experiment is now written to data.txt.
# - Writing of success.txt is more elegant now.
# - Implemented a new version of detectHbonds.
#
# Version 6.09:
# - Changed to a new version of mini-rsr that should improve pepflip's performance. This temporarily adds a lot of
#   dependencies.
# - The dependency problem will be solved once CCP4 starts distributing COOT 0.8.3.
# - Added a specific file for debug messages from the rebuilding tools in PDB_REDO.
#
# Version 6.08:
# - Failure of the second TLSANL run now triggers a debug message.
#
# Version 6.07:
# - Using a new version of pepflip that uses maps in which ligands an metals are masked out. This should improve the real-
#   space refinement in mini-rsr and reduce severe cases of severe backbone distortion.
# - Pepflip now uses the mini-rsr distributed with CCP4.
# - Fix to the second run of zen.
# - If different metal sites are discoverd in the second run of zen, the number of refinement cycles is increased.
# - Directories are now deleted without prompting.
#
# Version 6.06:
# - Now in k-fold cross validation, the coordinate perturbation is now very small and the number of refinement cycles is
#   increased.
# - Having more than one CRYST1 card now triggers a fatal error.
#
# Version 6.05:
# - SideAide can now use secondary structure specific rotamers. This is currently only done for Ile.
# - Stability fix for ion restraints.
# - Changed the tolerance for the density fit score in SideAide. This reduces the bias towards new conformations a bit.
# - Side chain flipping is now skipped for structures without protein. This gives a slight speedup for those structures.
# - Now properly deals with different tiers of YASARA.
# - UNL residues are now no longer deleted in local and in server mode.
#
# Version 6.04:
# - Extra metal, hydrogen bond, and nucleic acid restraints are re-evaluated after rebuilding.
#
# Version 6.03:
# - Jelly-body refinement can be switched off with the keyword '--nojelly'.
# - Removed the '-s' flag for picker in TLS group selection for better comparison of TLS group sections consisting of the
#   same number of groups.
# - Fixed the restraint pruning for stacking restraints.
#
# Version 6.02:
# - Zn-Cys4 clusters are now treated automatically. LINKs and angle restraints are added, disfulfide bridge detection in
#   REFMAC is switched off.
# - This option can be switched off with '--nometalrest'.
# - More verbose reporting on additional external restraints.
# - Hydrogen bond restraints can be used with '--hbondrest'.
#
# Version 6.01:
# - Unmerged symmetry related reflections are averaged using SFTOOLS rather than randomly selected by CAD.
# - The B-factor resetting now has a minumum of 10A^2.
# - The detection of disulfide bridges can be switched of with the hidden keyword --noss.
# - Moved the cispeptide detection to a separate YASARA run to improve stability.
# - Strict NCS and local NCS are no longer mutually exclusive.
# - Now doing a proper geometric weight optimisation in the lowest resolution category.
#
# Version 6.00:
# - First version with full OS X support.
# - Started reporting the solvent content from the density-based estimate from RWCONTENTS. Note that it is unreliable for
#   models with deuteriums, many alternates or strictncs.
# - Harmonic restraints are used for refinements with very small data sets (< 1000 reflections). They can be switched off
#   with the keyword '--noharmonic'.
# - Fix for generating dRSCC plots for really small structures.
# - Separated the restraints from LibG from user-supplied external restraints. The user-supplied restraints have a default
#   scale of 10, but this can be overruled in the restraint file using 'external weight scale [value]'
# - A warning is given when there are unmerged relections.
#
# Version 5.43:
# - Changes in the peptide isomer (cis to trans and vice versa) are listed in the COOT scripts.
# - Some minor syntax changes for the OSX version.
# - Improved LINK handeling.
#
# Version 5.42:
# - Switched to FoldX version 4.
# - Now writing compressed svg files.
# - A specific WHY_NOT message is written for PDB entries that only have an mmCIF coordinate file.
# - Rcomplete is also reported for full cross validation (courtesy of Tim Gruene).
# - Users can supply external restraints.
#
# Version 5.41:
# - Switched to CCP4-6.5.
# - LibG can be used to generate basepair and stacking restraints using the --libg command.
# - Removed TASER from the code.
# - Added a second pass of stripper after re-refinement to fix the LINK records.
# - In databank mode the FC, PHIC, FC_ALL_LS, and PHIC_ALL_LS are removed from the output mtz files.
#
# Version 5.40:
# - Better error message for models without CRYST1 records.
# - The killswitch is no longer used if the final refinement is restarted without updating the TLS.
# - Non-standard carbohydrate LINKs that involve the O1 atom are fixed by stripper.
# - Output from stripper is now written to the final log file.
#
# Version 5.39:
# - Bugfix in the detection of the original solvent model.
# - Chiron now fixes certain chirality problems by renaming a residue to its epimer.
# - Now used the default connectivity in settings Refmac (i.e. the 'connectivity' keyword was removed).
# - Changed the output webpage for the databank somewhat.
# - Only strict NCS is now used in rigid-body refinement.
# - Improved the detection of strict NCS.
# - SCALE cards with the identity matrix are now deleted.
#
# Version 5.38:
# - The boxplot and the dRSCC plot is now written as an svg file.
# - Moved the dRSCC plot generation to an external script.
# - B-factor model selection can also run in parallel now.
# - TLS group selection can also run in parallel now. The CPU limit to avoid hanging jobs now covers the entire selection
#   process and is increased to 24 hours.
# - Bond length and bond angle rmsZ scores that are so large that they cause character overflow are now also detected.
# - The real-space fit scores are now also stored for the initial PDB file.
#
# Version 5.37:
# - Stripper now deletes a number of LINK records that describe unlikely metal coordination.
# - Stripper now renames certain carbohydrates to improve the description of the biology.
# - Stripper writes Refmac style LINK records for certain carbohydrates to force the right LINK chirality.
# - Added multithreading for weight optimisations and cross validation. This gives a subsantial speed-up, but makes the
#   terminal output very ugly. The final log for server and local runs is cleaned up.
# - Multithreading required a rewrite of the error routines around Refmac. This makes them substantially lower and thus
#   PDB_REDO takes longer to stop with certain errors. Also problems with the TLS update after rebuilding are handled
#   differently (not much slower).
# - The second round of flipping is now only performed for structures with protein.
# - Removed some American spelling from the output.
#
# Version 5.36:
# - Stripper now rewrites LINKs for some carbohydrates.
# - Added the '--tighter' and '--looser' options to modify the restraint weight search space. This required rewriting the
#   way the restraints are set up.
#
# Version 5.35:
# - Compatibility fix for different versions of (g)awk.
# - The data.txt file now tells whether the data was treated as twinned.
# - Stability fix for the resolution cut-off determination.
# - Stability fix for reflection data conversion.
# - More verbose output for problems with extractor.
# - Implemented a fix for running PDB_REDO without YASARA.
# - B-factor resetting is now minimised at 5.00A^2.
# - Fixed residue selection for occupancy refinement.
#
# Version 5.34:
# - If the resolution cut-off is changed, RCAL and RFCAL are updated to represent the new data set.
# - Added a specific debug message for likely geometric restraint problems.
# - An update to binliner now devides extra reflection data into slightly more bins.
# - A number of fixes to extractor to deal with server input.
#
# Version 5.33:
# - Fixed a bug in the counting of residues with significant RSCC changes.
# - Fixed a bug in counting the number of peptide flips.
# - Added another specific error routine for poorly defined restraints.
# - Added two missing ILE rotamers to SideAide's search space.
# - In server and local mode, the log file of the last Refmac run is copied to the output.
# - If needed selenium atom occupancy in MSE is set to 1.00.
# - The success.txt file is now only written for the server.
# - Fixed a bug in generating the dRSCC plot for very small structures.
# - Added a reference to LigandExpo for misnamed atoms in local/server mode.
# - Added a warning message if cell dimensions from reflection data and model do not match in server/local mode.
# - Cif2cif now has an additional routine to check for suspicious sigF/sigI data.
# - Added the --maxres=VALUE keyword to cut the data at VALUE. The --cutreso flag was removed.
#
# Version 5.32:
# - Created a separate store for restraint files.
# - Fixed a bug for running on low-memmory mode.
# - Better handling of structures that need coordinate transformation.
# - Fixed a bug in setting the number of rigid body cycles.
# - The final model with total B-factors is now also written out.
#
# Version 5.31:
# - Added an error routine for when the FoldX licence has expired.
# - Resolute now compensates for missing free weighted R values.
# - Thanks to an update in freerflag twin-related reflections will now be part of the same test set.
# - Occupancy refinement is switched off for models that come straight from PHASER.
# - Explicit warning if kollumer cannot find usefull data columns.
#
# Version 5.30:
# - Several bugfixes in for TLS group extraction in extractor.
# - Slight changes for picker in extra strict mode.
# - Bugfix for extremely low resolution refinement.
# - Kollumer now also supports extracting anomalous intensities.
# - The identity MTRIX record is now added if needed in cases with strict NCS.
# - Improved error reporting, mostly for problems with the input data.
#
# Version 5.29:
# - Added a specific WHY_NOT message for insertion code specific NCS alignment problems.
# - Fixed a bug that caused unneeded validation of ions.
# - Added an error message for unmerged reflections in the input MTZ file.
#
# Version 5.28:
# - In the 'low' resolution group a tighter geometric restraint weight (.002) is also tested.
# - WHAT_TODO now also uses the output from extractor to prune the todo list.
# - There is now a second round of side-chain flipping if the side chains were rebuilt.
# - A guided tour through the model changes is written as a Coot script by a new program COOT_TOUR.
# - COOT_TOUR gives a more robust representation of H-bond flips, standardisation flips and changed rotamers.
# - Changed the rotamer definition so that residues with side-chain torsion angles only involving SP3 atoms require a
#   minimum torsion angle change of 60 degrees.
#
# Version 5.27:
# - PDB_REDO now makes a boxplot of the model quality with respect to the PDB and PDB_REDO.
# - Switched off the map generation. It is now handled by the databank server.
# - Added the '--notlsupdate' switch. This surpresses updating of the TLS tensors in the final refinement.
# - Bugfixes in generating the versions file.
# - Some more changes for OSX compatability.
# - If the I/sigI or F/sigF of the whole dataset is suspiciously low. It is assusmed that the sigF /sigI comlumn is fishy
#   and the column will be ignored. This required rewriting how cif2cif is run.
#
# Version 5.26:
# - PDB_REDO now writes a verbose summary in server mode.
# - Started using CCP4 6.4. This removes the dependency for DSSP and also requires support for the new output format from
#   EDSTATS.
# - Rewrote bits of code that used uniq, readlink or wget to add support for OSX.
# - Added back support for old WHAT_CHECK versions (only for OSX).
#
# Version 5.25:
# - Side chain rebuilding is now done for resoltions better than 3.3A. Peptide flipping from 3.5A and better. Other steps
#   are resolution independent.
# - Strict NCS matrices are now expanded during solvent mask parameter optimisation.
# - Stability fix for centrifuge.
#
# Version 5.24:
# - Stopped using experimental sigmas in the scaling function.
# - In server and local mode strict NCS restraints are switched on automatically.
#
# Version 5.23:
# - Stopped using the truncated mean for outlier rejection in the full cross-validation.
# - The descision to rebuild is now based on the numeric resolution (better than 3.5A) rather than on the resolution
#   category.
# - Improved the administration of whether or not a model was rebuilt. This solved a problem with the cross-validation.
#
# Version 5.22:
# - The TLS model generated by PDB_REDO will now always be called 'REDO.tls'.
# - Stopped making plots for dRSR. Adopted the output HTML accordingly.
# - Local and server runs now no longer give links to external websites.
# - HTML links to images and stylesheets now have absolute URLs to ensure fuctionality for local and server runs.
# - Only runs to the server will use the html version of the WHAT_CHECK reports. Local and databank runs will link to the
#   pdbout.txt files.
# - Other design updates to the HTML output.
# - A new WHAT_CHECK report will be generated for for each input model, even if a report is in the PDBREPORT databank.
#
# Version 5.21:
# - Atom coordinates are perturbed a bit before full cross validation if only one overall B-factor is used.
# - Improved reporting of use of NCS.
# - In local and server mode a Refmac command script is written.
#
# Version 5.20:
# - The --server flag now actually writes out specific files for the PDB_REDO server.
# - R-free is written out at more positions.
# - Output fix for the TLS group optimisation.
#
# Version 5.19:
# - The number of refinement cycles in the resolution test is the same as for other refinements (used to be always 30).
# - For every step of resolution increase the number of refinement cycles in in creased by five. New data is used so the
#   model should be further from convergence.
# - Added the '--nosigma' flag to not use the sigF or sigI column for scaling etc.
# - Do not try occupancy refinement for residues with an insertion code (not supported by Refmac).
#
# Version 5.18:
# - Modified amino acid residues are now also used at the termini of TLS groups.
# - Robustness update for the full cross validation.
# - A new file keeps track of unknown (computationally) chiral centers.
#
# Version 5.17:
# - Reimplemented the resolution cut-off optimisation with a tool called RESOLUTE.
# - When strictncs is used with a high number of copies (e.g. in viral capsids), the B-factor model is always isotropic
#   and reflections per atom are not used to set the resolution-based settings.
# - Separate message for problems parsing original TLS group selections.
#
# Version 5.16:
# - Changed definition of sigma(R-free): removed the factor two which inflates sigma(R-free), i.e. SIGRFCAL, SIGRFTLS, and
#   SIGRFFIN, by SQRT(2) and deflates the derived Z-scores by 1/SQRT(2).
# - Changed the cut-off for the R-free Z-score from 10 to 7 to detect possible R-free bias.
# - Changed the outlier rejection cut-off of longinus to Z = 2.6 (from 2.5).
# - Updated picker and changed its internal cut-off to 2.6 sigma (from 3.5).
# - Updated B-select as well.
#
# Version 5.15:
# - Lowered the cut-off for anisotropic B-factors to 30 reflections per atom.
# - Fixed debug message about missing sigma values.
# - Stability fix for mining WHAT_CHECK reports.
#
# Version 5.14:
# - Started using the new robust weighting scheme in Refmac.
# - Removed the '--noanisoscale' switch.
# - 2mFo-DFc maps in CCP4 format are generated for displaying purposes. Difference maps are genrated by the server.
# - No longer making a compressed file with the data. This is now done by the database's web server.
#
# A complete changelog is available from the PDB_REDO website
#
echo " "
if ($1 == "--local" || $1 == "--server") then
  setenv PDBID `perl -e '@chars=('a'..'z','A'..'Z','0'..'9'); $random_string; foreach (1..4) {$random_string.=$chars[rand @chars];} print "$random_string";'`
else
  setenv PDBID $1
  #Check the validity of the PDBid
  if ($PDBID != `echo $PDBID | cut -c 1-4`) then
    echo "This is not a valid PDB identifier. Cannot continue."
    exit(1)
  endif
endif
set D2       = `echo $PDBID | cut -c 2-3`    #Essential for certain environment values


# Check wether PDB_REDO is properly configured.
setenv TOOLS $BASE/tools                    #Directory with PDB_REDO software and data
if (! -e $TOOLS/pdb_redo.setup) then
  #PDB_REDO is not properly installed or configured.
  echo "FATAL ERROR!"
  echo "------------"
  echo "PDB_REDO cannot find its installation directory."
  echo "Please, set the correct value for 'BASE' in $0 on line 6."
  exit(1)
else
  #set additional environment parameters (define directories and files).
  source $TOOLS/pdb_redo.setup
endif

########################## Initialise important variables. Usually no need to edit this ##################################

#Add extra libraries
if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH $TOOLS/lib:$LD_LIBRARY_PATH
else
  setenv LD_LIBRARY_PATH $TOOLS/lib
endif

#Refinement parameters
set NCYCLE        = 20      #Initial number of (restrained) refinement cycles
set RBCYCLE       = 10      #Number of rigid-body refinement cycles
set BTESTCYC      = 10      #Number of restrained refinement cycles in the B-restraint weight optimisation
set TLSCYCLE      = 15      #Number of TLS-refinement cycles
set TLSCMD        =         #Empty TLS command
set TLSFILS       =         #No TLS files specified
set DOTLS         = 1       #Use TLS
set TLSUPDATE     = 1       #Update TLS tensors in final refinement
set OPTTLSG       = 'none'  #Fallback value if no TLS is used in refinement
set TIME          = `date +'%F'`
set BCMD          =         #Don't change the B-factors
set WGTTYPE       = MATRIX  #Use weight matrix by default
set BREFTYPE      = ISOT    #Do isotropic B-factor refinement by default
set NCSTYPE       =         #Empty NCS command
set RBNCS         =         #No NCS during rigid-body refinement
set NCSALIGN      =         #Only usefull when NCS is used (alignment cut-off)
set NCSNEIGH      =         #Only usefull when NCS is used (include neighbouring atoms in restraints)
set SOLVENT       = SIMP    #Use simple solvent model by default
set MASKPAR       =         #Use default masking parameters
set JELLY         =         #Do not use 'jelly body' refinement by default
set TORSION       =         #Do not use torsion restraints by default
set LOGSTEP       = 17      #Line number of Refmac log with required output (counted backwards)
set ESTRICT       =         #Picker is not extra strict
set BBONLY        = 0       #It is not a backbone only structure
set SHIFTCO       = 0.50    #Cut-off for atomic shift in ligand validation
set EXTSCAL       = 10      #Scale for external restraints
set REFIRES       =         #Resolution cut-offs for refinement; use all data by default
set OCCCMD        =         #Occupancy refinement commands
set METALCMD      =         #Command for metal site restraints
set RESTCMD       =         #Command for external restraints
set LIBGCMD       =         #DNA/RNA restraints commands
set HBONDWGTCMD   =         #Command for weighting external hydrogen bond restraints
set HBONDCMD      =         #Command for using hydrogen bond restraints
set HBRESTRWGT    = 2       #Hydrogen bond restraint weight
set NHBONDREST    = 0       #Number of H-bond restraints
set HOMOLWGTCMD   =         #Command for weighting external homology restraints
set HOMOLCMD      =         #Command for using homology restraints
set HOMOLRESTRWGT = 2       #Homology restraint weight
set NHOMOLREST    = 0       #Number of homology restraints
set KRABCMD       =         #Command for using restraints from KRAB
set KRABRESTRWGT  = 10      #Default restraint weight for KRAB restraints
set KRABWGTCMD    =         #Command for weighting KRAB restraints
set NKRABREST     = 0       #Number antibody-specific restraints 
set HODERCMD      =         #Command for extra homologous structures for hoder
set HARMCMD       =         #Harmonic restraints command
set WAVELCACHE    =         #Wavelength cashed from input MTZ file
set MAXPRUNE      = 20      #Give up after removing MAXPRUNE restraints
set SCATLIN       =         #Line that overruls the default scattering factor file
set SCATTERCMD    =         #Use the default X-ray scattering factors
set PHASES        =         #No phase columns as input by default
set ANOMCOEF      =         #Input anomalous data for Refmac
set ANOMCMD       =         #Command for anomalous refinement
set FASTAIN       =         #No input fasta by default

#Refmac settings
set CONNECTIVITY  = 'connectivity YES' #Empty value uses Refmac's default other options are 'connectivity [YES|NO|DEFINE]'
set SUGAR         =         #Empty value uses Refmac's default other options are 'sugar [YES|NO|DEFINE]'
set RSYMM         =         #Empty value uses Refmac's default other options are 'symmetry [YES|NO]'

#Default values
set NTLS         = 0       #Zero TLS groups by default. The number of groups is reset when TLS is used
set NPRUNEM      = 0       #No metal restraints pruned yet
set NPRUNEN      = 0       #No nucleic acid restraints pruned yet
set NPRUNE       = 0       #Total number of pruned restraints
set NMETALREST   = 0       #No metal restraints to start with
set NMETALREST2  = 0       #No metal restraints to start with
set EXPTYP       = 'X-ray' #The type of experiment
set GOTOLD       = 0       #By default not starting with previous PDB-REDO entry 

#Error flags
set TTEST      = 0 #No errors in the TLS group optimisation
set UTRUNCATE  = 0 #Truncate was not used

#Deal with potentially missing sigF/sigI values
set EXPSIG  = 'n'       #Do not use experimental sigmas for scaling
set WGTSIG  =           #Use experimental sigmas for scaling

#Debug flags and derivatives
set NOPDB        = 0
set NOSF         = 0
set USEMTZ       = 0
set DOWNLOAD     = 0
set NOHYD        = 0
set HYDROGEN     = ALL
set LEGACY       = 0
set DOTWIN       = 1
set DOHARMONIC   = 1      #Use harmonic restraints if needed
set TWIN         = 'test' #By default test for twinning
set ISTWIN       = 0      #Data not treated as twinned
set FALSETWIN    = 0      #Data not previously incrorrectly treated as twinned
set ISED         = 0      #this is an electron diffraction set
set DOTASER      = 0
set DOPEPFLIP    = 1
set DOSCBUILD    = 1
set DOCENTRIFUGE = 1
set DOREBUILD    = 1
set DONCS        = 1
set DORB         = 1 
set RESOCHECK    = 1
set DOJELLY      = 1     #Allow jelly-body refinement
set DOMETALREST  = 1     #Use special metal restraints by default
set STRICTNCS    = 0
set NCSSTRICT    =       #Empty strict NCS keyword for Refmac
set LOCAL        = 0
set SERVER       = 0
set NPROC        = 1     #Default number of processors to use
set MAXPROC      = 12    #The maximum number of processors to use
set PDBIN        =
set HKLIN        =
set MTZIN        =
set XTLS         =       #User submitted TLS files
set XHOM         =       #User submitted homologous structure models
set INREST       =       #The location of the geometric restraint file
set INEXT        =       #The location of the external restraint file
set INSEQ        =       #The location of the input sequence
set RELAX        =
set NEWMODEL     =
set LOWMEM       =
set FEWREFS      =
set CROSS        = 0
set C2CERR       = 0
set INTENS       =       #Use intensities if set to '-i'
set USTATUS      =       #Use the reflection status column. Set to '-s' to create a new R-free set
set SIGMA        =
set USIGMA       = 1     #Set to 0 to ignore sigma values
set C2CCONV      =       #Set to '-c' to convert intensities to amplitudes with cif2cif instead of ctruncate
set DOOCC        = 1     #Refine occupancies for certain residues
set DOLIBG       = 0
set CLEAN        = 0
set RESCATSTEP   = 0     #Do not change the resolution category
set SSBOND       = 'YES' #Autodetect disulfide bridges in Refmac
set HBONDREST    = 0
set DOKRAB       = 0     #Do not use antibody restraints by default
set DOHOMOLOGY   = 0     #Do not use homology restraints by default
set NOHOMOLOGY   = 0     #Modifier to switch off homology restraints
set SMODE        =       #Running mode for stripper
set C2CANO       = '-a'  #Run cif2cif in anomalous mode

#PDB-care settings
set CONPARAMS  =      #By default use the CONECT records

#Rebuilding settings
set NO_REBUILD =      #List of residues that must not be rebuilt
set H2O_KEEP   =      #List of waters to keep without question


############################################### You need not edit beyond this line #######################################

#OS specific settings:
if (`uname -s | grep -c arwin` != 0) then
  #OSX settings
  set ISOSX = 1
  alias readlink "python -c 'import os, sys; print os.path.realpath(sys.argv[2])'"
  setenv WEBGET 'curl -O -s'
else
  set ISOSX = 0
  setenv WEBGET 'wget -q'
  set    OLDWC = 0
endif

#Make sure the LOG file can be created (the :h modifier removes the filename from $LOG leaving only the directory)
mkdir -p $LOG:h

#make a temporary log file
if ($1 == "--local") then
  echo "Results from PDB_REDO $*" > $LOG
endif

#Write out header
echo " __   __   __      __   ___  __   __    __    _   _ " | tee -a $LOG
echo "|__) |  \ |__) __ |__) |__  |  \ /  \    /   / \ (_)" | tee -a $LOG
echo "|    |__/ |__)    |  \ |___ |__/ \__/   /  o \_/ (_)" | tee -a $LOG
echo " "

#Font for the header
# _     __  __       __  _  __  _   _
#/ \ /|  _) __) |_| |_  |_   / (_) (_|
#\_/  | /__ __)   | __) |_) /  (_)  _|


############################################### Check for debug flags ####################################################
#Get debug flags
echo " " | tee -a $LOG
echo "****** Optimisation setup ******" | tee -a $LOG
foreach ARG ($*)
  if ($ARG == $PDBID) then
    echo "-Evaluating PDB entry $PDBID" | tee -a $LOG
  else if ($ARG == "--local") then
    set LOCAL = 1
    #Set additional flags for programs
    set RELAX = "-r"  #Relaxed mode for extractor
    set SMODE = "-s"  #Server mode for stripper
  else if ($ARG == "--server") then
    set LOCAL  = 1
    set SERVER = 1
    set STDIR  = $cwd
    #Write different debug files
    set WHYNOT = $WHYNOT.server
    set DEBUG  = $DEBUG.server
    #Make the running status file
    touch $STDIR/processRunning.txt
    #Reset the logfile
    mkdir -p $STDIR/output
    echo "<pre style="\""font-family:'Lucida Console', Monaco, monospace"\"">" > $STDIR/output/process.log
    cat $LOG    >> $STDIR/output/process.log
    setenv LOG     $STDIR/output/process.log
    #Set additional flags for programs
    set RELAX = "-r"  #Relaxed mode for extractor
    set SMODE = "-s"  #Server mode for stripper
  else if (`echo $ARG | cut -c 1-8` == "--pdbin=") then
    #Get the file in two steps to expand relative paths, but also paths using '~'
    set PDBIN = `echo $ARG | cut -d '=' -f 2-`
    set PDBIN = `readlink -f $PDBIN`
    echo "-PDB_REDO will optimise the structure model in $PDBIN" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--hklin=") then
    set HKLIN = `echo $ARG | cut -d '=' -f 2-`
    set HKLIN = `readlink -f $HKLIN`
    #Choose the mmCIF file if two reflection files are specified
    if ($MTZIN != "") then
      echo "-You specified two reflection files. Using $HKLIN" | tee -a $LOG
      set MTZIN =
    else
      echo "-Using the experimental data in $HKLIN" | tee -a $LOG
    endif
  else if (`echo $ARG | cut -c 1-8` == "--mtzin=") then
    #Choose the mmCIF file if two reflection files are specified
    if ($HKLIN != "") then
      echo "-You specified two reflection files. Using $HKLIN" | tee -a $LOG
    else
      set MTZIN = `echo $ARG | cut -d '=' -f 2-`
      set MTZIN = `readlink -f $MTZIN`
      echo "-Using the experimental data in $MTZIN" | tee -a $LOG
    endif
  else if (`echo $ARG | cut -c 1-9` == "--dirout=") then
    set OUTPUT = `echo $ARG | cut -d '=' -f 2-`
    set OUTPUT = `readlink -f $OUTPUT`
  else if (`echo $ARG | cut -c 1-8` == "--tlsin=") then
    set INTLS = `echo $ARG | cut -d '=' -f 2-`
    set INTLS = `readlink -f $INTLS`
    echo "-Testing TLS group definitions from $INTLS" | tee -a $LOG
    #Push the input TLS file on a stack
    set XTLS  = `echo $XTLS $INTLS`
  else if (`echo $ARG | cut -c 1-8` == "--homin=") then
    set INHOM = `echo $ARG | cut -d '=' -f 2-`
    set INHOM = `readlink -f $INHOM`
    echo "-Using homologous structure model $INHOM" | tee -a $LOG
    #Push the input homologous file on a stack
    set XHOM  = `echo $XHOM $INHOM`    
  else if (`echo $ARG | cut -c 1-9` == "--restin=") then
    set INREST = `echo $ARG | cut -d '=' -f 2-`
    set INREST = `readlink -f $INREST`
    echo "-Using geometric restraints from $INREST" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--extin=") then
    set INEXT = `echo $ARG | cut -d '=' -f 2-`
    set INEXT = `readlink -f $INEXT`
    echo "-Using additional restraints from $INEXT" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--seqin=") then
    set INSEQ = `echo $ARG | cut -d '=' -f 2-`
    set INSEQ = `readlink -f $INSEQ`
    echo "-Using sequence from $INSEQ" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-8` == "--nproc=") then
    #Using multiple CPUs
    set NPROC = `echo $ARG | cut -d '=' -f 2-`
    #Is the value an integer?
    test $NPROC -eq $NPROC >& /dev/null
    if ($status) then
      #Not an integer, use a single core
      set NPROC = 1
      echo "-Using $NPROC CPU core" | tee -a $LOG
    else
      if ($NPROC < 1) then
        set NPROC = 1
      else if ($NPROC > $MAXPROC) then
        set NPROC = $MAXPROC
      endif
      echo "-Using a maximum of $NPROC CPU cores where possible" | tee -a $LOG
    endif
  else if ($ARG == --nopdb) then
    set NOPDB    = 1
    echo "-Keeping the PDB file in the working directory" | tee -a $LOG
  else if ($ARG == --nosf) then
    set NOSF     = 1
    echo "-Keeping the reflection file in the working directory" | tee -a $LOG
  else if ($ARG == --mtzdb) then
    set USEMTZ   = 1
    echo "-The reflection data in the data bank will be treated as MTZ" | tee -a $LOG
    set MTZIN = $SF/r${PDBID}sf.mtz
  else if ($ARG == --download) then
    set DOWNLOAD = 1
    echo "-The required files will be downloaded" | tee -a $LOG
  else if ($ARG == --nohyd) then
    set NOHYD    = 1
    set HYDROGEN = NO
    echo "-Hydrogens will NOT be added in the riding position during refinement" | tee -a $LOG
  else if ($ARG == --legacy) then
    set LEGACY = 1
    echo "-The PDB file will be treated as a legacy PDB entry" | tee -a $LOG
  else if ($ARG == --tighter) then
    @ RESCATSTEP = ($RESCATSTEP - 1)
    echo "-Trying tighter-than-default restraints" | tee -a $LOG
  else if ($ARG == --looser) then
    @ RESCATSTEP = ($RESCATSTEP + 1)
    echo "-Trying looser-than-default restraints" | tee -a $LOG
  else if ($ARG == --notls) then
    set DOTLS = 0
    echo "-No TLS refinement will be performed" | tee -a $LOG
  else if ($ARG == --notlsupdate) then
    set TLSUPDATE = 0
    echo "-No additional TLS refinement will be performed after rebuilding" | tee -a $LOG
  else if ($ARG == --noncs) then
    set DONCS = 0
    echo "-No NCS restraints will be used" | tee -a $LOG
  else if ($ARG == --nojelly) then
    set DOJELLY = 0
    echo "-Not doing jelly-body refinement" | tee -a $LOG
  else if ($ARG == --norb) then
    set DORB = 0
    echo "-Not doing rigid-body refinement" | tee -a $LOG    
  else if ($ARG == --noharmonic) then
    set DOHARMONIC = 0
    echo "-Not using harmonic restraints for small datasets" | tee -a $LOG
  else if ($ARG == --nometalrest) then
    set DOMETALREST = 0
    echo "-Not using additional metal site restraints" | tee -a $LOG
  else if ($ARG == --notwin) then
    set TWIN   =     #No detwinning
    set DOTWIN = 0
    echo "-No detwinning will be performed" | tee -a $LOG
  else if ($ARG == --noanomalous) then
    set C2CANO = ""
    echo "-Anomalous data will be ignored" | tee -a $LOG
  else if ($ARG == --newmodel) then
    set NEWMODEL = "-f" #Force rerefinement to finish with new model
    echo "-Rerefinement will always return a new model" | tee -a $LOG
  else if ($ARG == --nocentrifuge) then
    set DOCENTRIFUGE = 0    #No water deletion
    echo "-Poor waters will not be deleted" | tee -a $LOG
  else if ($ARG == --nopepflip) then
    set DOPEPFLIP = 0    #No peptide flips
    echo "-No peptide flips will be performed" | tee -a $LOG
  else if ($ARG == --noscbuild) then
    set DOSCBUILD = 0    #No side chain rebuilding
    echo "-Side chains will not be rebuilt" | tee -a $LOG
  else if ($ARG == --norebuild) then
    set DOREBUILD = 0    #No rebuilding at all
    echo "-No rebuilding at all" | tee -a $LOG
  else if ($ARG == --libg) then
    set DOLIBG = 1    #Generate basepair and stacking restraints
    echo "-Using DNA/RNA restraints from LibG" | tee -a $LOG
  else if ($ARG == --noocc) then
    set DOOCC = 0    #No occupancy refinement
    echo "-No occupancy refinement" | tee -a $LOG
  else if ($ARG == --notruncate) then
    set C2CCONV = '-c'    #Use cif2cif instead of ctruncate to convert intensities
    echo "-Not using truncate" | tee -a $LOG
  else if ($ARG == --lowmem) then
    #Reduce Refmac's memory consumption by not 'cleaning' the solvent mask
    set LOWMEM = "solvent process islands noremove"
    echo "-Working in limited memory mode" | tee -a $LOG
  else if ($ARG == --crossval) then
    set CROSS = 1
    echo "-Full k-fold cross validation will be performed" | tee -a $LOG
  else if ($ARG == --fewrefs) then
    set FEWREFS = "-n"
    echo "-Working with very small data set" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-9` == "--maxres=") then
    set MRESO = `echo $ARG | cut -d '=' -f 2-`
    echo "-The data will be cut to $MRESO A" | tee -a $LOG
  else if ($ARG == --intens) then
    set INTENS = "-i"
    echo "-Using intensities from the reflection file" | tee -a $LOG
  else if ($ARG == --nosigma) then
    set SIGMA  = "-g"
    set USIGMA = 0
    echo "-Not using sigF or sigI values from reflection file" | tee -a $LOG
  else if ($ARG == --noss) then
    set SSBOND = 'NO'
    echo "-Switched off detection of disulfide bridges" | tee -a $LOG
  else if ($ARG == --krab) then
    set DOKRAB = 1
    echo "-Using antibody-specific restraints" | tee -a $LOG    
  else if ($ARG == --hbondrest) then
    set HBONDREST = 1
    echo "-Using hydrogen bond restraints" | tee -a $LOG
  else if ($ARG == --homology) then
    set DOHOMOLOGY = 1
    echo "-Using homology-based restraints" | tee -a $LOG
  else if ($ARG == --nohomology) then
    set NOHOMOLOGY = 1
    echo "-No homology-based restraints will be used" | tee -a $LOG
  else if (`echo $ARG | cut -c 1-10` == "--comment=") then
    set COMMENT = `echo $* | grep -o -E -e '--comment='+'.*' | sed 's\--comment=\\' | sed 's\--.*\\'`
  else
    if (`echo $ARG | cut -c 1-2` == '--') then
      #It is an invalid keyword
      echo "Invalid keyword: $ARG" | tee -a $LOG
    else
      #Do nothing; it is part of the comment
    endif
  endif
end
if ($DOHOMOLOGY == 0 && "$XHOM" != "") then
  set DOHOMOLOGY = 1
  echo "-Using homology-based restraints" | tee -a $LOG
endif  

if ($LOCAL == 1) then
  echo "-During the process your structure model will be called '$PDBID'" | tee -a $LOG
  echo "-The final result will be written to $OUTPUT" | tee -a $LOG
endif

############################################### Parse the input files ####################################################

#For local runs make sure we have what we need
if ($LOCAL == 1) then
  #Is there a PDB file?
  if ($PDBIN == "") then
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "You did not provide an input PDB file. Cannot continue." | tee -a $LOG
    if ($SERVER == 0) then
      #Give a helpfull message
      echo "Use '--pdbin=some_pdb_file' and try again." | tee -a $LOG
    else
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
  #Is it a valid PDB file
  if (`grep -c CRYST1 $PDBIN` != 1) then
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "The input PDB file $PDBIN:t cannot be used. Please, provide a PDB file with one valid CRYST1 record." | tee -a $LOG
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
  #Is there reflection data
  if ($HKLIN == "" && $MTZIN == "") then
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    echo "You did not provide an input reflection file file. Cannot continue." | tee -a $LOG
    if ($SERVER == 0) then
      #Give a helpfull message
      echo "Use '--hklin=some_mmCIF_file' or '--mtzin=some_MTZ_file' and try again." | tee -a $LOG
    else
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
endif

#Delete previous results?
if ($NOPDB == 1 || $NOSF == 1) then
  cd $WORKDIR   #Go to workdir
  #Remove some files
  rm $WORKDIR/*.tls >& /dev/null
  rm -rf $WORKDIR/download >& /dev/null
  rm $WORKDIR/*.rest >& /dev/null
  rm $WORKDIR/mapval.log >& /dev/null
else
  #Delete all previous results!
  if (-e $WORKDIR) then
    rm -rf $WORKDIR
  endif
  mkdir -p $WORKDIR
  cd $WORKDIR
endif

#Report
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Parsing input files ******" | tee -a $LOG

#Download the files if the --download flag is given
if ($DOWNLOAD == 1) then
  echo "-Downloading" | tee -a $LOG
  #Make a download dir
  mkdir -p $WORKDIR/download
  cd $WORKDIR/download

  #Download the stuff (reflection data)
  #$WEBGET ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/structure_factors/r${PDBID}sf.ent.gz
  $WEBGET ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/structure_factors/r${PDBID}sf.ent.gz
  if($status) then
    echo " o Cannot download experimental data file" | tee -a $LOG
    exit(1)
  else
    echo " o Downloaded experimental data file" | tee -a $LOG
  endif
  gzip -df r${PDBID}sf.ent.gz
  setenv SF $WORKDIR/download

  #Download the stuff (model)
  #$WEBGET ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${PDBID}.ent.gz
  $WEBGET ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/all/pdb/pdb${PDBID}.ent.gz
  if($status) then
    echo " o Cannot download PDB file" | tee -a $LOG
    echo "COMMENT: No PDB-format coordinate file available" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
    exit(1)
  else
    echo " o Downloaded PDB file" | tee -a $LOG
  endif
  gzip -df pdb${PDBID}.ent.gz
  setenv PDB $WORKDIR/download

endif

#PDB file
if ($NOPDB == 0) then
  if ($LOCAL == 1) then
    if (-e $PDBIN) then
      #Copy PDB file
      cp $PDBIN $WORKDIR/pdb${PDBID}.ent
    else
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Cannot find the input file $PDBIN" | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  else
    #Copy PDB file from different directory structures
    if (-e $PDB/pdb${PDBID}.ent) then
      cp $PDB/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.ent
    else if (-e $PDB/$D2/pdb${PDBID}.ent) then
      cp $PDB/$D2/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.ent
    else if (-e $PDB/pdb${PDBID}.ent.gz) then
      zcat $PDB/pdb${PDBID}.ent.gz > $WORKDIR/pdb${PDBID}.ent      
    else if (-e $PDB/$D2/pdb${PDBID}.ent.gz) then
      zcat $PDB/$D2/pdb${PDBID}.ent.gz > $WORKDIR/pdb${PDBID}.ent
    else
      echo "-No PDB-format coordinate file." | tee -a $LOG
      #Test for the existence of a reflection file
      echo "COMMENT: No PDB-format coordinate file available" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
      cd $BASE
      exit(1)
    endif
  endif
endif

#Does the input PDB file have SEQRES records?
if (`grep -c ^SEQRES $WORKDIR/pdb${PDBID}.ent` > 0) then
  set GOTSEQRES = 1
else
  set GOTSEQRES = 0
endif

#Check the experimental method (first through the EXPDTA record then by looking for other clues)
if (`grep -c EXPDTA $WORKDIR/pdb${PDBID}.ent` == 0) then
  set ISXRAY = 1 #Assume X-ray and do not check experimental method further
else
  set ISXRAY = `grep EXPDTA $WORKDIR/pdb${PDBID}.ent | head -n 1 | grep -c X-RAY`
endif
if ($ISXRAY == 0) then #Not an X-ray structure, is it electron diffraction?
  set ISED = `grep EXPDTA $WORKDIR/pdb${PDBID}.ent | head -n 1 | grep -c 'ELECTRON CRYSTALLOGRAPHY'`
  if ($ISED == 1) then
    echo "-Switching to electron diffraction mode" | tee -a $LOG
    #Use converted X-ray to electron form factors
    set SCATTERCMD = 'source EM MB'
    set EXPTYP     = 'ED'
  else
    #Create a WHY_NOT comment
    echo "-This is not an X-ray diffraction structure" | tee -a $LOG
    echo "COMMENT: Not an X-ray structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#If this is a multi-entry complex, do not re-refine, create WHY_NOT comment
if (`grep -c -E '^SPLIT' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This is part of a structure divided over multiple PDB entries" | tee -a $LOG
  echo "COMMENT: Multi-entry complex" >> $WHYNOT
  echo "PDB-REDO,$PDBID"              >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Check for implicit strict NCS
if (`grep -E '^MTRIX' $WORKDIR/pdb${PDBID}.ent | grep -c -v -E '^MTRIX'+'.{54}1'` != 0) then
  set STRICTNCS = 1
  set NCSSTRICT = ncsconstraints
endif

#If implicit atoms are used as described in REMARK 285, create WHY_NOT comment
if (`grep -c -E '^REMARK 285 THE ENTRY PRESENTED HERE DOES NOT CONTAIN THE COMPLETE' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This PDB file has implicit atoms" | tee -a $LOG
  echo "COMMENT: Asymmetric unit not complete" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                       >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#If the atoms are not in the 'standard crystal frame' according to REMARK 285, give a warning
if (`grep -c -E 'ARE NOT PRESENTED IN THE STANDARD CRYSTAL' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-The atoms in this PDB are not in the standard crystal frame" | tee -a $LOG
  #Write a warning iff a transformation is required
  if (`grep -c 'IN ORDER TO GENERATE THE CRYSTAL AU' $WORKDIR/pdb${PDBID}.ent` != 0) then
    echo " o Coordinate transformation needed. Cannot continue." | tee -a $LOG
    echo "COMMENT: Atoms not in standard crystal frame" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                              >> $WHYNOT
    exit(1)
  endif
endif

#If it is a multi-model refinement structure, create WHY_NOT comment
if (`grep -c -E '^ENDMDL' $WORKDIR/pdb${PDBID}.ent` != 0) then
  echo "-This is a multi-model structure" | tee -a $LOG
  echo "COMMENT: Multi-model refinement structure" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                           >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Check to see if it is a CA-only structure
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` == 1) then
  #There is only one type of atom. check to see if it is CA
  if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "CA"` != 0) then
    echo "-This is a C-alpha-only structure" | tee -a $LOG
    echo "COMMENT: C-alpha-only structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Check to see if it is a CA/P-only structure
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` == 2) then
  #There is only two types of atom. check to see if they are CA and P
  if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "CA"` != 0 && \
      `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "P "` != 0) then
    echo "-This is a C-alpha/backbone-phosphorus-only structure" | tee -a $LOG
    echo "COMMENT: C-alpha/backbone-phosphorus-only structure" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                     >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Check to see if it is a backbone-only structure
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` == 4) then
  #There are only four types of atom. Check to see if they are the backbone atoms.
  if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "CA"` != 0 && \
      `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "C "` != 0 && \
      `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "N "` != 0 && \
      `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | grep -c "O "` != 0) then
    set BBONLY = 1
    #Check for peptide-like hetero compounds with side chains (residues with more than 4 atoms)
    foreach RES (`grep -E '^HETATM'+'.{6}'+'( CA| C | N | O )' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | sort -u`)
      if (`grep -E '^HETATM'+'.{11}'+$RES $WORKDIR/pdb${PDBID}.ent | cut -c 13-16 | sort -u | wc -l` > 4) then
        #This is not a backbone-only structure
        set BBONLY = 0
      endif
    end
    if($BBONLY == 1)then
      echo "-This is a backbone-only structure" | tee -a $LOG
      echo "COMMENT: Backbone-only structure" >> $DEBUG
      echo "PDB-REDO,$PDBID"                  >> $DEBUG
    endif
  endif
endif

#Is there only 1 type of residue and is it UNK?
if (`grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | sort -u | wc -l` == 1 && \
    `grep -E '^ATOM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | grep -c "UNK"` != 0) then
  #Do not use autoNCS (poly-UNKs cannot be alligned properly)
  echo "-The main chain consists only of UNK residues, NCS will not be used!" | tee -a $LOG
  set DONCS = 0
  #Are there ligands or hetero compounds (count the number of different residues)
  if (`grep -E '^ATOM|^HETATM' $WORKDIR/pdb${PDBID}.ent | cut -c 18-20 | sort -u | wc -l` == 1) then
    echo "-This is an UNK-only structure" | tee -a $LOG
    echo "COMMENT: UNK-only structure" >> $DEBUG
    echo "PDB-REDO,$PDBID"             >> $DEBUG
  endif
endif

#Is it a PDB file that came straight from PHASER?
if (`grep -c 'REMARK Log-Likelihood Gain' $WORKDIR/pdb${PDBID}.ent` != 0 && `grep -c 'BUSTER' $WORKDIR/pdb${PDBID}.ent` == 0) then
  #Switch off occupancy refinement
  set DOOCC = 0

  #Report and write a debug message
  echo "-PHASER output model detected. Switching off occupancy refinement." | tee -a $LOG
  echo "COMMENT: PHASER model; no occupancy refinement" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                >> $DEBUG
endif

#Get temporary cell diemensions and space group
set MAAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $2}'`
set MBAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $3}'`
set MCAXIS = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $4}'`
set MALPHA = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $5}'`
set MBETA  = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $6}'`
set MGAMMA = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | awk '{print $7}'`
set MSPACE = `grep ^CRYST1 $WORKDIR/pdb${PDBID}.ent | cut -c 56-66 | sed "s/P 21 21 2 A/P 21 21 2 (a)/g" | sed "s/P 1-       /P -1/g"`


#Structure factors
if ($NOSF == 0) then
  #Get the user supplied reflection file...
  if ($LOCAL == 1 || $USEMTZ == 1) then
    #If supplied, analyse the input reflection file and convert it to mmCIF if needed...
    if ($MTZIN != "") then
      if (-e $MTZIN) then
        #Dump the MTZ file. PROGRAM: mtzdump
        mtzdmp $MTZIN > $WORKDIR/cifcreate.log
        if ($status) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-Cannot parse the input MTZ file." | tee -a $LOG
            echo "COMMENT: Cannot parse the input MTZ file" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                          >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "The input reflection file $MTZIN:t is not in the MTZ format." | tee -a $LOG
            echo "Please, rerun PDB_REDO with reflections in the MTZ format." | tee -a $LOG
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Stop if we are dealing with an unmerged reflection file
        if (`grep -c 'Y  M/ISYM' $WORKDIR/cifcreate.log` != 0) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-The input MTZ file has unmerged reflections. Cannot continue." | tee -a $LOG
            echo "COMMENT: Input MTZ file has unmerged reflections" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "The MTZ file $MTZIN:t seems to contain unmerged reflections." | tee -a $LOG
            echo "Please, rerun PDB_REDO with merged reflections." | tee -a $LOG
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Cache the wavelength
        set WAVELCACHE = `grep -A 6 'Dataset ID, project/crystal/dataset names, cell' $WORKDIR/cifcreate.log | tail -n 1`

        #Sanity check for the wavelength cache
        if ($WAVELCACHE == 0.00000) then
          set WAVELCACHE =       #Empty value
        endif

        #Get column labels PROGRAM: kollumer
        set LABELS = `$TOOLS/kollumer $WORKDIR/cifcreate.log`

        #Warn if no proper columns were found
        if (`echo $LABELS | grep -c Problem` != 0) then
          #Give error message
          if ($USEMTZ == 1) then
            #Local databank-type message
            echo "-No useful information extracted from MTZ file. Cannot continue." | tee -a $LOG
            echo "COMMENT: No useful data in input MTZ file" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                           >> $WHYNOT
          else
            #--local type error
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Could not extract useful reflection data from $MTZIN:t" | tee -a $LOG
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
          endif
          cd $BASE
          exit(1)
        endif

        #Convert the MTZ file to mmCIF. PROGRAM: mtz2various
        mtz2various \
        HKLIN  $MTZIN \
        HKLOUT $WORKDIR/r${PDBID}sf.ent \
        <<eof >> $WORKDIR/cifcreate.log
          OUTPUT CIF data_$PDBID
          LABIN $LABELS
          END
eof

        #Read the cell dimensions for validation
        set RAAXIS = `grep _cell.length_a $WORKDIR/r${PDBID}sf.ent | awk '{printf "%.3f\n", $2}'`
        set RBAXIS = `grep _cell.length_b $WORKDIR/r${PDBID}sf.ent | awk '{printf "%.3f\n", $2}'`
        set RCAXIS = `grep _cell.length_c $WORKDIR/r${PDBID}sf.ent | awk '{printf "%.3f\n", $2}'`
        set RALPHA = `grep _cell.angle_alpha $WORKDIR/r${PDBID}sf.ent | awk '{printf "%.2f\n", $2}'`
        set RBETA  = `grep _cell.angle_beta $WORKDIR/r${PDBID}sf.ent  | awk '{printf "%.2f\n", $2}'`
        set RGAMMA = `grep _cell.angle_gamma $WORKDIR/r${PDBID}sf.ent | awk '{printf "%.2f\n", $2}'`

        #Give a warning if the cell dimensions do not match
        if ($RAAXIS != $MAAXIS || $RBAXIS != $MBAXIS || $RCAXIS != $MCAXIS || $RALPHA != $MALPHA || $RBETA != $MBETA || $RGAMMA != $MGAMMA) then
          #The warning
          echo " " | tee -a $LOG
          echo "WARNING!" | tee -a $LOG
          echo "--------" | tee -a $LOG
          echo "The cell dimensions from the reflection data and the model do not match!" | tee -a $LOG
          echo "From reflections: $RAAXIS $RBAXIS $RCAXIS $RALPHA $RBETA $RGAMMA" | tee -a $LOG
          echo "From model      : $MAAXIS $MBAXIS $MCAXIS $MALPHA $MBETA $MGAMMA" | tee -a $LOG
          echo "Cell dimensions from the input model will be used." | tee -a $LOG
          echo " " | tee -a $LOG
          #Write debug message
          echo "COMMENT: cell dimension mismatch" >> $DEBUG
          echo "PDB-REDO,$PDBID"                  >> $DEBUG
        endif

      else
        #Give error message (--local type)
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot find the reflection file $MTZIN" | tee -a $LOG
        if ($SERVER == 1) then
          #Write out status files
         touch $STDIR/stoppingProcess.txt
         touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    else
      #...or just copy the mmCIF file
      if (-e $HKLIN) then
        cp $HKLIN $WORKDIR/r${PDBID}sf.ent
      else
        #Give error message (--local type)
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot find the reflection file $HKLIN" | tee -a $LOG
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    endif
  #...or the reflection file from (a local copy of) the PDB
  else
    if (-e $SF/r${PDBID}sf.ent) then
      cp $SF/r${PDBID}sf.ent $WORKDIR/
    else if (-e $SF/$D2/r${PDBID}sf.ent.gz) then
      cp $SF/$D2/r${PDBID}sf.ent.gz $WORKDIR/
      gzip -df $WORKDIR/r${PDBID}sf.ent.gz
    else if (-e $SF/r${PDBID}sf.ent.gz) then
      cp $SF/r${PDBID}sf.ent.gz $WORKDIR/
      gzip -df $WORKDIR/r${PDBID}sf.ent.gz
    else
      echo "-No structure factors found" | tee -a $LOG
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  endif
else
  #Clean some stuff up to avoid problems later on
  if (-e $WORKDIR/$PDBID.cif) then
    rm $WORKDIR/$PDBID.cif
  endif
endif


#Create WHAT_CHECK report
echo "-Validating the input PDB file with WHAT_CHECK" | tee -a $LOG

#Check the WHAT_CHECK version (if it is below 10, we need different settings)
if ($ISOSX == 1) then
  $WC/bin/whatcheck \
  << eof >& $WORKDIR/versions.log
eof
  set OLDWC = `grep version $WORKDIR/versions.log | cut -c 34-37 | awk '{if ($1 < 10) {print "1"} else {print "0"}}'`
endif

#Initialise
set WCERR = 0

#Go to temporary running directory
setenv WCWORO $WORKDIR/wctemo
mkdir -p $WCWORO
cd $WCWORO

#Get the PDB file
cp $WORKDIR/pdb$PDBID.ent $WCWORO

#Do the actual validation. PROGRAM: WHAT_CHECK
$WC/bin/whatcheck $WCWORO/pdb$PDBID.ent Y >& $WORKDIR/wc_original.log

#Check for an output file
if (-e $WCWORO/pdbout.txt) then
  #Do nothing
else
  #Give warning
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate original structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                        >> $DEBUG
endif

#Create webpage
$WC/bin/pdbout2html >>& $WORKDIR/wc_original.log

#Check index.html completeness
if (`grep -c "Final summary" $WCWORO/pdbout.html` == 0 && `grep -c "Summary report" $WCWORO/pdbout.html` == 0) then
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate original structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                        >> $DEBUG
endif

#Clean up on aisle six
mkdir -p $WORKDIR/wo
set PDBOUT = $WORKDIR/wo/pdbout.txt

if (-e $WCWORO/pdbout.txt) then
  mv $WCWORO/pdbout.txt $WORKDIR/wo/
endif
if (-e $WCWORO/check.db) then
  bzip2 -f $WCWORO/check.db
  mv $WCWORO/check.db.bz2 $WORKDIR/wo/
endif
#Copy html and images only for the server
if ($SERVER == 1) then
  mv $WCWORO/pdbout.html $WORKDIR/wo/index.html
  mv $WCWORO/*.gif $WORKDIR/wo/ >& /dev/null
endif

#Go to start directory
cd $WORKDIR
rm -rf $WCWORO


################################################# Remove unwanted atoms ##################################################

echo "-Preparing the structure model" | tee -a $LOG

renumbered:

#Perform DEFY flips before reannotating any LINKs. PROGRAM: flipper
echo " o Performing DEFY flips" | tee -a $LOG
cp $WORKDIR/pdb$PDBID.ent $WORKDIR/pdb$PDBID.bak
$TOOLS/flipper -v \
-pdbin  $WORKDIR/pdb$PDBID.bak \
-pdbout $WORKDIR/pdb$PDBID.ent > $WORKDIR/flipper.log

renumbered:

#Run PDB-care if it is defined
if ($?PDBCARE) then
  echo " o Running pdb-care" | tee -a $LOG
  #Modify behaviour based on the presence/absence of CONECT records
  if (`grep -c '^CONECT' $WORKDIR/pdb$PDBID.ent` == 0) then
    set CONPARAMS = "-p $UMF_HOME/conectparams.xml"
  endif

  #Run pdb-care. PROGRAM: umfconverter
  $PDBCARE/umfconverter \
  -i pdb $WORKDIR/pdb$PDBID.ent \
  -o $WORKDIR/$PDBID.val \
  -f c \
  $CONPARAMS \
  > $WORKDIR/pdbcare.log
  if ($status) then
    echo "   * Error using pdb-care" | tee -a $LOG
    echo "COMMENT: umfconverter: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                      >> $DEBUG
  endif
endif

#Check whether PDB-care file exists
if (-e $WORKDIR/$PDBID.val) then
  set CAREXML = $WORKDIR/$PDBID.val
else
  set CAREXML =        #No PDB-care XML file
endif

#Check for the presence of Zn atoms if allowed
if ($DOMETALREST == 1 && `grep '^[AH][TE][OT][MA]' $WORKDIR/pdb$PDBID.ent | cut -c 18-20 | grep -E -c 'ZN | ZN'` > 0) then
  
  #Make temporary mtz file
  unique \
  HKLOUT $WORKDIR/temp.mtz \
<<eof > $WORKDIR/mtz_creationt.log
  CELL $MAAXIS $MBAXIS $MCAXIS $MALPHA $MBETA $MGAMMA
  SYMM '$MSPACE'
  LABOUT F=FUNI SIGF=SIGFUNI
  RESOLUTION 12.0
  END
eof

  #Run Zen. PROGRAM: zen
  echo " o Running Zen" | tee -a $LOG
  cp $WORKDIR/pdb$PDBID.ent $WORKDIR/pdb$PDBID.older
  $TOOLS/zen \
  -v $WORKDIR/pdb$PDBID.older \
  $WORKDIR/pdb$PDBID.ent \
  $WORKDIR/metal.raw \
  $WORKDIR/temp.mtz \
  > zen.log

  #Set up the use of the restraints
  if (`cat $WORKDIR/metal.raw | wc -l` > 0) then
    set METALCMD = "@$WORKDIR/metal.rest"
  endif

  #Cleanup
  sed 's/ \+/ /g' $WORKDIR/metal.raw >> $WORKDIR/metal.rest
  rm $WORKDIR/metal.raw
  rm $WORKDIR/temp.mtz
endif

#Delete atoms
echo " o Running stripper" | tee -a $LOG

#Make backup
cp $WORKDIR/pdb$PDBID.ent $WORKDIR/pdb$PDBID.old

#Run stripper to remove all hydrogens, deuteriums, atoms with type X, unknown ligands (UNL), UNK atoms not described in
#the refmac library, and superfluous carbohydrate oxygen atoms. Also remove crazy LINKs. PROGRAM: stripper
$TOOLS/stripper -v $SMODE \
$WORKDIR/pdb$PDBID.old \
$WORKDIR/pdb$PDBID.ent \
$TOOLS/pdb_redo.dat \
$CAREXML \
> $WORKDIR/stripper.log

#Report changes
echo "   * `tail -n 6 $WORKDIR/stripper.log | head -n 1`" | tee -a $LOG
echo "   * `tail -n 5 $WORKDIR/stripper.log | head -n 1`" | tee -a $LOG
echo "   * `tail -n 4 $WORKDIR/stripper.log | head -n 1`" | tee -a $LOG
echo "   * `tail -n 3 $WORKDIR/stripper.log | head -n 1`" | tee -a $LOG
echo "   * `tail -n 2 $WORKDIR/stripper.log | head -n 1`" | tee -a $LOG
echo "   * `tail -n 1 $WORKDIR/stripper.log`" | tee -a $LOG


#Correct some things and run stripper again
if (`tail -n 1 $WORKDIR/stripper.log | awk '{print $4}'` != 0 || `grep 'Added   LINK   records' $WORKDIR/stripper.log | awk '{print $4}'` != 0) then
  #Use pdbcur to sort the atoms if carbohydrate LINKs were fixed. PROGRAM: pdbcur
  if (`tail -n 1 $WORKDIR/stripper.log | awk '{print $4}'` != 0) then
    cp $WORKDIR/pdb$PDBID.ent $WORKDIR/pdb$PDBID.cur

    pdbcur \
    XYZIN $WORKDIR/pdb$PDBID.cur \
    XYZOUT $WORKDIR/pdb$PDBID.ent \
    <<eof > $WORKDIR/pdbcur.log
    END
eof
  endif

  #Run PDB-care again if needed
  if (`grep 'Added   LINK   records' $WORKDIR/stripper.log | awk '{print $4}'` != 0) then
    #Store old .val file
    mv $WORKDIR/$PDBID.val $WORKDIR/$PDBID.val1

    #Run PDB-care
    $PDBCARE/umfconverter \
    -i pdb $WORKDIR/pdb$PDBID.ent \
    -o $WORKDIR/$PDBID.val \
    -f c \
    $CONPARAMS \
    >> $WORKDIR/pdbcare.log
    if ($status) then
      echo "   * Error using pdb-care" | tee -a $LOG
      echo "COMMENT: umfconverter: general error 2" >> $DEBUG
      echo "PDB-REDO,$PDBID"                        >> $DEBUG
    endif
  endif

  #Store old pdb file
  cp $WORKDIR/pdb$PDBID.ent $WORKDIR/pdb$PDBID.old2

  #Run stripper again
  $TOOLS/stripper -v $SMODE \
  $WORKDIR/pdb$PDBID.old2 \
  $WORKDIR/pdb$PDBID.ent \
  $TOOLS/pdb_redo.dat \
  $CAREXML \
  >> $WORKDIR/stripper.log
endif

########################## Extract essential information from structure factors and pdb file #############################

#Check structure factor file and reformat
#Label to come back to for forced use of intensities
beintens:

#Run cif2cif
echo "-Checking reflection data" | tee -a $LOG
c2cagain:

#PROGRAM: cif2cif
$TOOLS/cif2cif $C2CANO $USTATUS $FEWREFS $INTENS $SIGMA $C2CCONV \
$WORKDIR/r${PDBID}sf.ent \
$WORKDIR/$PDBID.cif \
$WAVELCACHE \
> $WORKDIR/${PDBID}c2c.log

#Success or not?
if (-e $WORKDIR/$PDBID.cif) then

  #Check for status flag column (only once)
  if ($?GOTR) then
    #GOTR is already set. Do nothing.
  else
    #GOTR is 1 if the status flag (R/Rfree) is used, 0 if not
    set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
    if ($GOTR == 0) then
      #Check the reason for the missing set
      if (`grep -c 'Warning: too small R-free set.' $WORKDIR/${PDBID}c2c.log` != 0) then
        echo " o Original R-free set too small, a new one will be created" | tee -a $LOG
      else
        echo " o No R-free set defined, a new one will be created" | tee -a $LOG
      endif
    endif
  endif
else if (! $?GOTR) then
  #...warn about the R-free set and try again
  echo " o Error(s) in structure factors; not using status flag" | tee -a $LOG
  set C2CERR = 1
  set GOTR = 0
  set USTATUS = '-s'
  mv $WORKDIR/${PDBID}c2c.log $WORKDIR/${PDBID}c2cv1.log
  goto c2cagain
else
  #Something is really wrong. Halt and give WHY_NOT comment
  echo " o Cannot parse structure factor file. Exit." | tee -a $LOG
  echo "COMMENT: cif2cif: cannot parse structure factors" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
  rm core.*
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  exit(1)
endif

#Check for fishy sigF or sigI data
if (`grep -c 'using the -g switch!' $WORKDIR/${PDBID}c2c.log` != 0) then
  #Give a warning
  echo " o Suspicious sigma values detected; they will be ignored" | tee -a $LOG
  echo "COMMENT: cif2cif: suspicious sigma values" >> $DEBUG
  echo "PDB-REDO,$PDBID"                           >> $DEBUG

  #Go back to cif2cif, now ignorig the sigma column.
  mv $WORKDIR/${PDBID}c2c.log $WORKDIR/${PDBID}c2cv2.log
  mv $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_badsig.cif
  set SIGMA = '-g'
  goto c2cagain
endif

#Check for multiple wavelength data
if (`grep -c "Second wavelength dataset found!" $WORKDIR/${PDBID}c2c.log` != 0) then
  echo " o Only using reflection data from the first wavelength" | tee -a $LOG
  echo "COMMENT: cif2cif: only using data from first wavelength" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                         >> $DEBUG
endif

#Check to see whether a suitable set of experimental sigmas exists (checked now because DYEAR is needed)
if (`grep -c "sigma values" $WORKDIR/${PDBID}c2c.log` != 0) then
  set EXPSIG  =  'n' #Do not use experimental sigmas for weighting
  set WGTSIG  = "NOEX"
  set USIGMA  = 0
  if (`grep -c "No experimental sigmas found" $WORKDIR/${PDBID}c2c.log` != 0) then
    echo " o No usable experimental sigmas in reflection data"    | tee -a $LOG
    #Give debug message later
  else
    echo " o Cannot use experimental sigmas for weighting" | tee -a $LOG
    #Give debug message later
  endif
endif

#Check for phase information (only for electron diffraction)
if ($ISED == 1) then
  if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0) then
    echo " o Experimental phases in HL format will be used"     | tee -a $LOG
    set PHASES = "HLA=HLA HLB=HLB HLC=HLC HLD=HLD"
    set TWIN   =     #No detwinning
    set DOTWIN = 0
  else if (`grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
    echo " o Experimental phases with figures of merit will be used" | tee -a $LOG
    set PHASES = "PHIB=PHIB FOM=FOM"
    set TWIN   =     #No detwinning
    set DOTWIN = 0
  endif
endif

#Report on anomalous data
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  if (`echo $PHASES | cut -c 1` == "") then
    echo " o Anomalous data will be used to calculate anomalous maps" | tee -a $LOG
  endif
endif

#Extract information from a pdbfile en the newly created reflection file
#TLS group information is extracted from the PDB header. If no groups are defined, every chain gets its own TLS group
echo "-Extracting data from the PDB file" | tee -a $LOG

#PROGRAM: extractor
$TOOLS/extractor -f $RELAX \
$WORKDIR/pdb${PDBID}.ent \
$WORKDIR/${PDBID}.cif \
$CLIBD_MON/list/mon_lib_list.cif \
$TOOLS/pdb_redo.dat \
$WORKDIR/$PDBID.extracted \
$WORKDIR/$PDBID.tls > $WORKDIR/extractor.log

#Success or not?
if (-e $WORKDIR/$PDBID.extracted) then
#Do nothing
else
  #Give PDB-REDO mode-specific error messages
  if ($LOCAL == 1) then
    #Give the long error message
    echo " " | tee -a $LOG
    echo " " | tee -a $LOG
    echo "FATAL ERROR!" | tee -a $LOG
    echo "------------" | tee -a $LOG
    if (`grep -c 'No space group in the PDB file' $WORKDIR/extractor.log` != 0) then
      echo "Space group missing in the PDB file. Please add it." | tee -a $LOG
    else if (`grep -c 'Cannot interpret the PDB file:' $WORKDIR/extractor.log` != 0) then
      echo "Cannot intrepret the input PDB file at this line:" | tee -a $LOG
      grep -A 1 'Cannot interpret the PDB file:' $WORKDIR/extractor.log | tail -n 1 | tee -a $LOG
      echo "Please, ensure that you provide a valid PDB file."
    else
      echo "Could not parse the input PDB file. Please, ensure it is a valid file." | tee -a $LOG
    endif
  else
    #Give the short error message
    echo " " | tee -a $LOG
    echo " o Cannot parse pdb file. Exit." | tee -a $LOG
  endif
  echo "COMMENT: extractor: error using PDB file" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                          >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Are there any TLS groups from extractor?
if (-e $WORKDIR/$PDBID.tls || -e $WORKDIR/REDO.tls) then
  if (`grep -c 'TLS groups created' $WORKDIR/extractor.log` != 0) then
    echo " o `grep 'TLS groups created' $WORKDIR/extractor.log`" | tee -a $LOG
  endif
  if (`grep -c 'TLS groups extracted' $WORKDIR/extractor.log` != 0) then
    echo " o `grep 'TLS groups extracted' $WORKDIR/extractor.log`" | tee -a $LOG
  endif
else
  echo " o Extractor could not create TLS group definitions" | tee -a $LOG
  echo "COMMENT: extractor: cannot create TLS groups" >> $DEBUG
  echo "PDB-REDO,$PDBID"                              >> $DEBUG
endif

#Check for TLS parsing errors
if (`grep -c "negative selectors" $WORKDIR/extractor.log` != 0 || `grep -c "Cannot read TLS residue range" $WORKDIR/extractor.log` != 0) then
  echo " o Extractor could not parse original TLS group definitions" | tee -a $LOG
  echo "COMMENT: extractor: cannot parse TLS groups" >> $DEBUG
  echo "PDB-REDO,$PDBID"                             >> $DEBUG
endif

#Copy in extra TLS files
if ("$XTLS" != "") then
  echo "-Importing extra TLS group definitions" | tee -a $LOG
  set CNT = 0
  #Loop over files
  foreach FIL ($XTLS)
    set CNT  = `expr $CNT + 1`
    set RANK = `seq -w $CNT 99 | head -n 1`
    #Only copy the group definitions, not the tensors
    grep -E 'TLS|RANGE|^ ' $FIL > $WORKDIR/in$RANK.tls
    #Reject unusable files
    if (`grep -c 'RANGE' $WORKDIR/in$RANK.tls` == 0) then
      echo " o The file $FIL does not contain valid TLS group definitions" | tee -a $LOG
      rm $WORKDIR/in$RANK.tls
      set CNT  = `expr $CNT - 1`
    endif
  end
  echo " o Definitions imported: $CNT" | tee -a $LOG
endif

#Get important parameters from $PDBID/$PDBID.extracted
set AAXIS      = `head -n 1  $WORKDIR/$PDBID.extracted`
set BAXIS      = `head -n 2  $WORKDIR/$PDBID.extracted | tail -n 1`
set CAXIS      = `head -n 3  $WORKDIR/$PDBID.extracted | tail -n 1`
set ALPHA      = `head -n 4  $WORKDIR/$PDBID.extracted | tail -n 1`
set BETA       = `head -n 5  $WORKDIR/$PDBID.extracted | tail -n 1`
set GAMMA      = `head -n 6  $WORKDIR/$PDBID.extracted | tail -n 1`
set RESOLUTION = `head -n 7  $WORKDIR/$PDBID.extracted | tail -n 1`
setenv RESO $RESOLUTION  #Set value for bits of Perl code in the script
set DATARESH   = `head -n 8  $WORKDIR/$PDBID.extracted | tail -n 1`
setenv DRESH $DATARESH   #Set value for bits of Perl code in the script
set DATARESL   = `head -n 9  $WORKDIR/$PDBID.extracted | tail -n 1`
set RFACT      = `head -n 10 $WORKDIR/$PDBID.extracted | tail -n 1`
setenv PRFACT $RFACT
set RFREE      = `head -n 11 $WORKDIR/$PDBID.extracted | tail -n 1`
setenv PRFREE $RFREE
set NO_REBUILD = `head -n 12 $WORKDIR/$PDBID.extracted | tail -n 1`
set BAVER      = `head -n 13 $WORKDIR/$PDBID.extracted | tail -n 1`
setenv PBAVER $BAVER
#set BREF       = `head -n 14 $WORKDIR/$PDBID.extracted | tail -n 1` #Not used
set REFCNT     = `head -n 15 $WORKDIR/$PDBID.extracted | tail -n 1`
set TSTCNT     = `head -n 16 $WORKDIR/$PDBID.extracted | tail -n 1`
set TSTPRC     = `head -n 17 $WORKDIR/$PDBID.extracted | tail -n 1`
set PROG       = `head -n 18 $WORKDIR/$PDBID.extracted | tail -n 1`
set DYEAR      = `head -n 19 $WORKDIR/$PDBID.extracted | tail -n 1`
set SPACEGROUP = `head -n 20 $WORKDIR/$PDBID.extracted | tail -n 1`
set H2O_KEEP   = `head -n 21 $WORKDIR/$PDBID.extracted | tail -n 1`
set BBN_KEEP   = `head -n 22 $WORKDIR/$PDBID.extracted | tail -n 1`
set BBO_KEEP   = `head -n 23 $WORKDIR/$PDBID.extracted | tail -n 1`
set GOT_PROT   = `head -n 24 $WORKDIR/$PDBID.extracted | tail -n 1`
set VDWPROBE   = `head -n 25 $WORKDIR/$PDBID.extracted | tail -n 1`
set IONPROBE   = `head -n 26 $WORKDIR/$PDBID.extracted | tail -n 1`
set RSHRINK    = `head -n 27 $WORKDIR/$PDBID.extracted | tail -n 1`
set COMPLETEH  = `head -n 28 $WORKDIR/$PDBID.extracted | tail -n 1`
set LIG_LIST   = `head -n 29 $WORKDIR/$PDBID.extracted | tail -n 1`
set GOT_NUC    = `head -n 30 $WORKDIR/$PDBID.extracted | tail -n 1`
set WAVELPDB   = `tail -n 1 $WORKDIR/$PDBID.extracted`
setenv PCOMPLETEH $COMPLETEH

#Check for strict NCS and count the number of atoms
#Calculate the number of atoms in the refinement...
setenv PATMCNT `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/pdb${PDBID}.ent`
setenv PMTRIX "1"
#... if strict NCS is used, multiply the apparent number of atoms with the number of MTRIX records
if ($STRICTNCS == 1) then
  #Count the number of MTRIX records
  setenv PMTRIX `grep -E '^MTRIX' $WORKDIR/pdb${PDBID}.ent | cut -c 7-10 | sort -u | wc -l`
  #Get the total
  setenv PATMCNT `perl -e '$totatom = $ENV{PMTRIX}*$ENV{PATMCNT}; print "$totatom";'`
endif

#Get the solvent content
rwcontents \
XYZIN $WORKDIR/pdb${PDBID}.ent <<eof > $WORKDIR/rwcont.log
eof
if ($status) then
  #Give a debug message
  echo "COMMENT: rwcontents failed"   >> $DEBUG
  echo "PDB-REDO,$PDBID"              >> $DEBUG
  set SOLVD = 'NA'
else
  #Get the numbers
  set SOLVD = `grep 'Assuming protein density is' $WORKDIR/rwcont.log | cut -d ':' -f 2`
endif

#Get the sequence
if ($INSEQ != "") then
  echo " o Importing amino acid sequence" | tee -a $LOG
  cp $INSEQ $WORKDIR/user.fasta
  set FASTAIN = "-fastain $WORKDIR/user.fasta"
endif

#Extract and or check the sequence
if ($GOT_PROT == 'T') then
  echo " o Extracting the amino acid sequence" | tee -a $LOG
  $TOOLS/pdb2fasta \
  -v \
  $FASTAIN \
  -pdb $WORKDIR/pdb${PDBID}.ent \
  -fasta $WORKDIR/$PDBID.fasta \
  -tools $TOOLS \
  >& $WORKDIR/pdb2fasta.log
  #make sure there is a fasta file
  if (! -e $WORKDIR/$PDBID.fasta) then
    if ($DOHOMOLOGY == 1) then
      echo "  * Cannot extract sequence. Not using homology restraints." | tee -a $LOG
      set DOHOMOLOGY = 0
    else
      echo "   * Cannot extract sequence." | tee -a $LOG
    endif
    echo "COMMENT: Cannot make fasta file" >> $DEBUG
    echo "PDB-REDO,$PDBID"                 >> $DEBUG
  else if (-e $WORKDIR/user.fasta) then
    #There is a fasta file check if there were sequence conflicts
    if (`grep -c 'mis-matched' $WORKDIR/pdb2fasta.log` > 0) then
      #Show the sequince conflicts
      echo " " | tee -a $LOG
      echo "WARNING!" | tee -a $LOG
      echo "--------" | tee -a $LOG
      echo "Conflict(s) between the input sequence and the input model detected." | tee -a $LOG
      echo "Please, check whether these conflicts can be resolved." | tee -a $LOG
      echo " " | tee -a $LOG
      echo "Details from pdb2fasta output:" | tee -a $LOG
      grep 'mis-matched' $WORKDIR/pdb2fasta.log | cut -c 10- | tee -a $LOG
      echo "--------" | tee -a $LOG
      echo " " | tee -a $LOG
    endif
  endif
endif

#Check to see whether a suitable set of experimental sigmas exists (checked now because DYEAR is needed)
#Only give error message is the structure was deposited this century
if ($DYEAR > 2009 && `grep -c "sigma values" $WORKDIR/${PDBID}c2c.log` != 0 ) then
  if (`grep -c "No experimental sigmas found" $WORKDIR/${PDBID}c2c.log` != 0 && $USIGMA == 1) then
    echo "COMMENT: cif2cif: no expertimental sigmas found"   >> $DEBUG
    echo "PDB-REDO,$PDBID"                                   >> $DEBUG
  endif
endif

#Debug message for missing R-free set
if ($C2CERR == 1 && $DYEAR > 2009 && `grep -c useful $WORKDIR/${PDBID}c2c.log` == 0) then
  if (`grep -c useful $WORKDIR/${PDBID}c2c.log` == 0) then
    echo "COMMENT: cif2cif: cannot use _refln.status column" >> $DEBUG
    echo "PDB-REDO,$PDBID"                                   >> $DEBUG
  endif
endif

#Set mask parameter
set MASKPAR  = "solvent vdwprobe $VDWPROBE ionprobe $IONPROBE rshrink $RSHRINK"

#Check for freakishly high resolution
if (`perl -e 'if ($ENV{DRESH} < 0.30) {print "1"} else {print "0"};'` == 1) then
  echo " o Unlikely high resolution. Reflection data may be corrupted. Cannot continue." | tee -a $LOG
  echo "COMMENT: Suspiciously high data resolution" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                            >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

#Set resolution cut-offs if needed
if (`perl -e 'if ($ENV{RESO} - $ENV{DRESH} > 0.10) {print "1"} else {print "0"};'` == 1) then
  set REFIRES = "RESO $RESO"
  setenv URESO $RESO
else
  setenv URESO $DATARESH
endif

#Check for missing high resolution data
if (`perl -e 'if ($ENV{DRESH} - $ENV{RESO} > 0.10) {print "1"} else {print "0"};'` == 1) then
  if ($LOCAL == 1) then
    echo "WARNING!" | tee -a $LOG
    echo "--------" | tee -a $LOG
    echo "The resolution of the model is higher than that of the data." | tee -a $LOG
    echo "Please, check whether the correct reflection data file was used." | tee -a $LOG
  else
    echo " o High resolution data is missing, calculated R-factors are unrealiable" | tee -a $LOG
  endif
  #Write debug warning
  echo "COMMENT: High resolution data missing" >> $DEBUG
  echo "PDB-REDO,$PDBID"                       >> $DEBUG
  
  #Compensate
  if ($?COMMENT) then
    set COMMENT = "$COMMENT; High resolution data is missing" 
  else
    set COMMENT = "High resolution data is missing" 
  endif
  set NEWMODEL = "-f" #Force rerefinement to finish with new model
endif

#Set legacy mode for PDB entries from the seventies and eighties and ED structures predating 1995
if ($DYEAR < 1990 && $ISXRAY == 1) then
  set LEGACY = 1
  echo " o Warning pre-1990 X-ray PDB entry. It will be treated as legacy entry." | tee -a $LOG
  echo "COMMENT: Run in legacy mode" >> $DEBUG
  echo "PDB-REDO,$PDBID"             >> $DEBUG
else if ($DYEAR < 1995 && $ISED == 1) then
  set LEGACY = 1
  echo " o Warning pre-1995 electron diffraction PDB entry. It will be treated as legacy entry." | tee -a $LOG
  echo "COMMENT: Run in legacy mode" >> $DEBUG
  echo "PDB-REDO,$PDBID"             >> $DEBUG
endif

#Stop and give WHY_NOT comment if no R-factor can be extracted from the PDB header...
if ($RFACT == 0.9990) then
  # ... unless PDB-REDO is running in legacy mode or working on a local file
  if ($LEGACY == 1 || $LOCAL == 1) then
    echo " o Cannot extract R-factor from PDB header. Recalculated value will be used as refinement target." | tee -a $LOG
    if ($LOCAL == 1) then
      set LEGACY = 1
      echo " o Switching to legacy mode." | tee -a $LOG
    endif
  else
    echo " o Cannot extract R-factor from PDB header. Cannot continue." | tee -a $LOG
    echo "COMMENT: No R-factor reported" >> $WHYNOT
    echo "PDB-REDO,$PDBID"               >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Obtain solvent model from REMARK records
if (`grep '^REMARK   3' $WORKDIR/pdb${PDBID}.ent | grep 'METHOD USED' | grep -c -E 'BABINET|SWAT|BULK|MOEWS|KRETSINGER|TNT'` != 0) then
  set SOLVENT = BULK
endif
#Fall back to a simple model when the keywords 'FLAT' or 'CNS BULK' are found.
if (`grep '^REMARK   3' $WORKDIR/pdb${PDBID}.ent | grep 'METHOD USED' | grep -c -E 'MASK|FLAT|CNS BULK'` != 0) then
  set SOLVENT = SIMP
endif

#Strip out the PDBID form the header record
set UPID = `echo $PDBID | awk '{print toupper($0)}'`
cp  $WORKDIR/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.tmp
sed "s/ $UPID / XXXX /g" $WORKDIR/pdb${PDBID}.tmp > $WORKDIR/pdb${PDBID}.ent
rm  $WORKDIR/pdb${PDBID}.tmp

############################# Create MTZ file for structure factor handling in Refmac ####################################

echo "-Creating MTZ file" | tee -a $LOG

mtzmaking:

#Import CIF file (using all reflections). PROGRAM: cif2mtz
cif2mtz \
HKLIN  $WORKDIR/$PDBID.cif \
HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
  CELL  $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA
  SYMM '$SPACEGROUP'
  END
eof
if ($status) then
  echo " o Error using CIF2MTZ. Cannot continue." | tee -a $LOG
  echo "COMMENT: cif2mtz: general error" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                 >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  exit(1)
endif

#Remove the phase columns if there is anomalous data and we do not need to do phased refinement
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  #We have anomalous data, do we also have phases?
  if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
    #Do we need phased refinement?
    if (`echo $PHASES | cut -c 1` == "") then
      #Make a back-up
      cp $WORKDIR/raw.mtz $WORKDIR/raw_phase.mtz

      #Which labels must be removed
      if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0 && `grep -c '_refln.phase_meas' $WORKDIR/$PDBID.cif` != 0) then
        set DELLABEL = "HLA HLB HLC HLD PHIB FOM"
      else if (`grep -c '_refln.pdbx_HL_A_iso' $WORKDIR/$PDBID.cif` != 0) then
        set DELLABEL = "HLA HLB HLC HLD"
      else
        set DELLABEL = "PHIB FOM"
      endif

      #Strip out the phase columns
      mtzutils \
      HKLIN  $WORKDIR/raw_phase.mtz \
      HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
        EXCLUDE $DELLABEL
        END
eof

    endif
  endif
endif

#Cut the resolution if needed
if ($?MRESO) then
  echo " o Cutting the data at $MRESO A" | tee -a $LOG
  set DATARESH   = $MRESO
  setenv DRESH $DATARESH   #Set value for bits of Perl code in the script
  setenv URESO $DATARESH

  cp $WORKDIR/raw.mtz $WORKDIR/raw_uncut.mtz
  #Run MTZUTILS to cut the data. PROGRAM: mtzutils
  mtzutils \
  HKLIN  $WORKDIR/raw_uncut.mtz \
  HKLOUT $WORKDIR/raw.mtz \
  <<eof >>$WORKDIR/mtz_creation.log
    RESOLUTION $DATARESH $DATARESL
    END
eof
endif

#Add the wavelength if possible
#Returnpoint
wavelengthadd:

if (`grep -c '_diffrn_radiation_wavelength.wavelength' $WORKDIR/$PDBID.cif` != 0) then
  #Get the wavelenghth
  set WAVELENGTH = `grep '_diffrn_radiation_wavelength.wavelength' $WORKDIR/$PDBID.cif | awk '{print $2}'`
else if ($WAVELPDB != '0.00000') then
  set WAVELENGTH = $WAVELPDB
else
  set WAVELENGTH = 'NA'
  if ($DYEAR > 2006) then
    echo "COMMENT: wavelength unknown" >> $DEBUG
    echo "PDB-REDO,$PDBID"             >> $DEBUG
  endif
endif

#Add the wavelength to the MTZ file
if ($WAVELENGTH != 'NA') then
  #Make a back-up
  cp $WORKDIR/raw.mtz $WORKDIR/raw_nowavel.mtz

  #Run CAD. PROGRAM: cad
  cad \
  HKLIN1 $WORKDIR/raw_nowavel.mtz \
  HKLOUT $WORKDIR/raw.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 1  ALLIN
    DWAVELENGTH FILE_NUMBER 1 1 $WAVELENGTH
    END
eof
  #Check for a large reduction in the number of reflections. This indicates unmerged reflections.
  setenv PREFIN  `grep 'and a total of' $WORKDIR/mtz_creation.log | tail -n 1 | awk '{print $5}'`
  setenv PREFOUT `grep 'Final Total of Unique records to HKLOUT' $WORKDIR/mtz_creation.log | tail -n 1 | cut -d '=' -f 2`
  if (`perl -e 'if ($ENV{PREFIN} / $ENV{PREFOUT} > 1.9) {print "1"} else {print "0"};'` == 1) then
    #There seem to be unmerged reflections or Friedel pairs on separate lines or many systematically absent reflections
    if (`grep -c 'Systematic absent reflection rejected' $WORKDIR/mtz_creation.log` > 1000) then
      #Something seroiusly wrong with the dataset
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "There are too many systematic absent reflections. Cannot continue. Please, check the input reflection data." | tee -a $LOG
      else
        echo " o There are too many systematic absent reflections. Cannot continue." | tee -a $LOG
      endif
      echo "COMMENT: Too many systematic absent reflections" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      exit(1)
    endif

    #Report only once
    if ($?MERGEREP) then
      #Problem already reported
    else
      echo " o The reflections are not completely merged" | tee -a $LOG
      echo "COMMENT: unmerged reflections" >> $DEBUG
      echo "PDB-REDO,$PDBID"               >> $DEBUG
      set MERGEREP = 1
    endif

    #First get the sigF or sigI column if it was rejected before
    if (-e $WORKDIR/${PDBID}_badsig.cif) then
      #Report
      echo "   * Recovering sigma value data before merging" | tee -a $LOG

      #Take the reflection file with the sigma columns
      mv $WORKDIR/${PDBID}_badsig.cif $WORKDIR/$PDBID.cif

      #Remake the mtz file
      goto mtzmaking
    endif

    #Delete the output from cad.
    rm $WORKDIR/raw.mtz

    #Merge the reflections. PROGRAM: sftools
    sftools << eof >> $WORKDIR/mtz_creation.log
    read $WORKDIR/raw_nowavel.mtz
    reduce
    merge average
    write $WORKDIR/raw.mtz
    stop
eof

    #Add the wavelength again
    goto wavelengthadd
  endif
endif

#create a unique list of reflections for the given unit cell-symmetry-resolution. PROGRAM: unique
unique \
HKLOUT $WORKDIR/unique.mtz \
<<eof >> $WORKDIR/mtz_creation.log
  CELL $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA
  SYMM '$SPACEGROUP'
  LABOUT F=FUNI SIGF=SIGFUNI
  RESOLUTION $DATARESH
  END
eof

#Combine data (generate R-free set if one doesn't already exist)
if ($GOTR == 0) then

  #Calculate the required R-free fraction (between 0.05 and 0.10)
  setenv PREFCNT $REFCNT  #Set value for bits of Perl code in the script
  set FRAC = `perl -e 'if ($ENV{PREFCNT} > 20000) {print "0.05"} elsif ($ENV{PREFCNT} < 10000) {print "0.10"} else {printf ("%.4f\n", 1000/$ENV{PREFCNT})};'`

  #Generate R-free flag, then ....
  echo " o Generating R-free set using $FRAC fraction of all reflections" | tee -a $LOG

  #PROGRAM: freerflag
  freerflag \
  HKLIN $WORKDIR/unique.mtz \
  HKLOUT $WORKDIR/freer.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    FREERFRAC $FRAC
    END
eof

  #.... combine data
  cad \
  HKLIN2 $WORKDIR/freer.mtz \
  HKLIN1 $WORKDIR/raw.mtz \
  HKLOUT $WORKDIR/combined.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 2  E1=FreeR_flag
    LABIN FILE 1  ALLIN
    END
eof

  #Rename FreeR_flag to FREE
  mtzutils \
  HKLIN  $WORKDIR/combined.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    COLUMN_LABELS FreeR_flag=FREE
    END
eof

  #Clean up extra temporary files
  rm $WORKDIR/freer.mtz

else

  #Just combine data
  cad \
  HKLIN2 $WORKDIR/unique.mtz \
  HKLIN1 $WORKDIR/raw.mtz \
  HKLOUT $WORKDIR/combined.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    LABIN FILE 2  ALLIN
    LABIN FILE 1  ALLIN
    END
eof

  #Strip out the FUNI and SIGFUNI columns
  mtzutils \
  HKLIN  $WORKDIR/combined.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    EXCLUDE FUNI SIGFUNI
    END
eof
endif

#Clean up
rm $WORKDIR/combined.mtz

#Need to work with intensities?
if (`mtzdmp $WORKDIR/merged.mtz -e | grep -A 2 'Column Types' | grep -c J` != 0) then
  echo " o Using ctruncate to convert intensities to amplitudes" | tee -a $LOG
  set UTRUNCATE = 1
  #Create a back-up
  cp $WORKDIR/merged.mtz $WORKDIR/mergedbu.mtz

  if (`grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
    #Convert I to F, also for anomalous data. PROGRAM: ctruncate
    ctruncate \
    -mtzin  $WORKDIR/merged.mtz \
    -mtzout $WORKDIR/ctruncate.mtz \
    -no-aniso \
    -freein "/*/*/[FREE]" \
    -colano "/*/*/[I(+),SIGI(+),I(-),SIGI(-)]" \
    -colin  "/*/*/[I,SIGI]" >>& $WORKDIR/mtz_creation.log
  else
    #Convert I to F. PROGRAM: ctruncate
    ctruncate \
    -mtzin  $WORKDIR/merged.mtz \
    -mtzout $WORKDIR/ctruncate.mtz \
    -no-aniso \
    -freein "/*/*/[FREE]" \
    -colin "/*/*/[I,SIGI]" >>& $WORKDIR/mtz_creation.log
  endif

  if ($status || ! -e $WORKDIR/ctruncate.mtz) then
    #Try running without using ctruncate
    echo "   * Error using ctruncate" | tee -a $LOG
    echo "   * Using cif2cif to convert intensities" | tee -a $LOG

    #Make backup
    cp $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_notruncate.cif

    #Run cif2cif
    $TOOLS/cif2cif -c $C2CANO $FEWREFS $INTENS $USTATUS $SIGMA \
    $WORKDIR/${PDBID}_notruncate.cif \
    $WORKDIR/$PDBID.cif \
    >> $WORKDIR/mtz_creation.log

    #Check for the existence of a test set.    
    set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
    if ($GOTR == 0) then
      #Check the reason for the missing set
      if (`grep -c 'Warning: too small R-free set.' $WORKDIR/mtz_creation.log` != 0) then
        echo "   * Original R-free set too small, a new one will be created" | tee -a $LOG
      else
        echo "   * No R-free set defined, a new one will be created" | tee -a $LOG
      endif
    endif
    
    #Restart making the MTZ file
    goto mtzmaking
  endif

  #Rename amplitude columns and remove intensity columns
  if (`grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
    set DELLABEL = 'DANO SIGDANO ISYM I(+) SIGI(+) I(-) SIGI(-)'
  else
    set DELLABEL = 'I SIGI'
  endif

  mtzutils \
  HKLIN  $WORKDIR/ctruncate.mtz \
  HKLOUT $WORKDIR/merged.mtz \
<<eof >> $WORKDIR/mtz_creation.log
    COLUMN_LABELS F=FP SIGF=SIGFP
    EXCLUDE $DELLABEL
    END
eof

  #Clean up
  rm $WORKDIR/ctruncate.mtz
endif


#Clean up data and create output file (use separate log file)
freerflag \
HKLIN $WORKDIR/merged.mtz \
HKLOUT $WORKDIR/$PDBID.mtz \
<<eof >& $WORKDIR/freerflag.log
  COMPLETE FREE=FREE
  END
eof
if ($status) then
  #Problem with FREERFLAG
  echo " o Problem with test set, creating a new one." | tee -a $LOG

  #Give debug message
  echo "COMMENT: freerflag: test set problem" >> $DEBUG
  echo "PDB-REDO,$PDBID"                      >> $DEBUG

  #Strip old test set
  mtzutils \
  HKLIN  $WORKDIR/merged.mtz \
  HKLOUT $WORKDIR/nofree.mtz \
<<eof >> $WORKDIR/freerflag.log
    EXCLUDE FREE
    END
eof

  #Calculate the required R-free fraction (between 0.05 and 0.10)
  setenv PREFCNT $REFCNT  #Set value for bits of Perl code in the script
  set FRAC = `perl -e 'if ($ENV{PREFCNT} > 20000) {print "0.05"} elsif ($ENV{PREFCNT} < 10000) {print "0.10"} else {printf ("%.4f\n", 1000/$ENV{PREFCNT})};'`

  #Generate R-free flag, then ....
  echo "   * Generating R-free set using $FRAC fraction of all reflections" | tee -a $LOG

  freerflag \
  HKLIN  $WORKDIR/nofree.mtz \
  HKLOUT $WORKDIR/newfree.mtz \
<<eof >> $WORKDIR/freerflag.log
    FREERFRAC $FRAC
    END
eof

  #Rename FreeR_flag to FREE
  mtzutils \
  HKLIN  $WORKDIR/newfree.mtz \
  HKLOUT $WORKDIR/$PDBID.mtz \
<<eof >> $WORKDIR/freerflag.log
    COLUMN_LABELS FreeR_flag=FREE
    END
eof

endif

#Append freerflag.log to mtz_creation.log
cat $WORKDIR/freerflag.log >> $WORKDIR/mtz_creation.log
rm  $WORKDIR/freerflag.log


#Create lowres mtz file for sfcheck and mtz2various
if (`perl -e 'if ($ENV{RESO} - $ENV{DRESH} > 0.10) {print "1"} else {print "0"};'` == 1) then
  #Run MTZUTILS to cut the data
  mtzutils \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/lowres.mtz \
  <<eof >> $WORKDIR/mtz_creation.log
    RESOLUTION $RESOLUTION $DATARESL
    END
eof
  if ($status) then
    echo "   * Error using MTZUTILS. Cannot continue." | tee -a $LOG
    echo "COMMENT: mtzutils: general error" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                  >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif
else
  cp $WORKDIR/$PDBID.mtz $WORKDIR/lowres.mtz
endif

#Calculate Wilson B, completeness and twinning stats. PROGRAM: SFcheck
sfcheck \
-f $WORKDIR/lowres.mtz \
-mem 150 \
>> $WORKDIR/mtz_creation.log

#See if sfcheck ran propely and truncate was used
if (`grep -c "ERR: NUMBER OF REFLS" sfcheck.log` != 0 && $UTRUNCATE == 1) then
  #Try running without using ctruncate
  echo "   * Error using ctruncate" | tee -a $LOG
  echo "   * Using cif2cif to convert intensities" | tee -a $LOG

  #Make backup
  cp $WORKDIR/$PDBID.cif $WORKDIR/${PDBID}_notruncate.cif

  #Run cif2cif
  $TOOLS/cif2cif -c $C2CANO $FEWREFS $INTENS $USTATUS $SIGMA \
  $WORKDIR/${PDBID}_notruncate.cif \
  $WORKDIR/$PDBID.cif \
  >> $WORKDIR/mtz_creation.log
  
  #Check for the existence of a test set.    
  set GOTR = `grep -c _refln.status $WORKDIR/$PDBID.cif`
  if ($GOTR == 0) then
    #Check the reason for the missing set
    if (`grep -c 'Warning: too small R-free set.' $WORKDIR/mtz_creation.log` != 0) then
      echo "   * Original R-free set too small, a new one will be created" | tee -a $LOG
    else
      echo "   * No R-free set defined, a new one will be created" | tee -a $LOG
    endif
  endif

  #Restart making the MTZ file
  goto mtzmaking
endif

#Mine sfcheck output
setenv BWILS     `grep 'by Wilson' sfcheck.log | awk '{print $7}'`
setenv COMPLETED `grep 'Completeness :' sfcheck.log | cut -c 16- | awk '{print $1}'`
setenv PTWINA    `grep 'Alpha(twin fraction)' sfcheck.log | cut -c 38-43`

#Check for twinning with Phaser if SFCHECK suggests twinning. PROGRAM: PHASER
if ($PTWINA != "") then
  phaser \
  << eof > $WORKDIR/phaser.log
    MODE NCS
    ROOT $WORKDIR
    HKLIN $WORKDIR/lowres.mtz
    LABIN  F=FP SIGF=SIGFP
    COMPOSITION BY SOLVENT
    COMPOSITION PERCENTAGE $SOLVD
    HKLOUT OFF
    END
eof
  if (`grep -c "EXIT STATUS: SUCCESS" $WORKDIR/phaser.log` == 0) then
    #Phaser problem.
    echo "COMMENT: phaser: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                >> $DEBUG
    set PHASTWIN = 0
  else
    #Check whether PHASER found twinning
    if (`grep -c "Warning: Intensity moments suggest significant twinning" $WORKDIR/phaser.log` == 0 && `grep -c "Warning: Intensity moments suggest possibility of twinning" $WORKDIR/phaser.log` == 0) then
      #Phaser did not find twinning
      set PHASTWIN = 0
    else
      echo " o Both SFCHECK and PHASER suggest the data are twinned" | tee -a $LOG
      set PHASTWIN = 1
    endif
  endif
else
  #Don't look for twinning in Phaser
  set PHASTWIN = 0
endif


#Get the correct number of (R-free) reflections
mtz2various \
HKLIN $WORKDIR/lowres.mtz \
HKLOUT $WORKDIR/temp.hkl \
<<eof >> $WORKDIR/mtz_creation.log
  OUTPUT CIF data_temp
  LABIN FP=FP SIGFP=SIGFP FREE=FREE
  END
eof

set NTSTCNT = `grep -cE ' f ' $WORKDIR/temp.hkl`
setenv PNTSTCNT $NTSTCNT  #Set value for bits of Perl code in the script
set NREFCNT = `grep -cE ' [of] ' $WORKDIR/temp.hkl`
setenv PWORKCNT `grep -cE ' o ' $WORKDIR/temp.hkl`

#Check for anomalous data
if (`grep -c '_refln.pdbx_F_plus' $WORKDIR/$PDBID.cif` != 0 || `grep -c '_refln.pdbx_I_plus' $WORKDIR/$PDBID.cif` != 0) then
  set ANOMCOEF = 'F+=F(+) SIGF+=SIGF(+) F-=F(-) SIGF-=SIGF(-)'
  set ANOMCMD  = 'ANOM maponly'
endif


#Cleanup
rm sfcheck.xml
rm sfcheck_XXXX.ps
rm $WORKDIR/temp.hkl
rm $WORKDIR/lowres.mtz


###################################### Report structure model and data details ###########################################

#Print values to screen. Compensate if R(-free) was not reported.
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure model and data details ******" | tee -a $LOG
echo "Cell axes            : $AAXIS $BAXIS $CAXIS" | tee -a $LOG
echo "Cell angles          : $ALPHA $BETA $GAMMA"  | tee -a $LOG
echo "Space group          : $SPACEGROUP" | tee -a $LOG
echo "Resolution (pdb)     : $RESOLUTION" | tee -a $LOG
echo "Resolution (data)    : $DATARESH"   | tee -a $LOG  #Data resolution (only reflections with F/SIGF>1.0 are used)
echo "Lowest resolution    : $DATARESL"   | tee -a $LOG #Lowest resolution in the dataset (only reflections with F/SIGF>1.0 are used)
if ($RFACT == 0.9990) then
  echo "Reported Rfactor     : none"   | tee -a $LOG
  set RFACT = "NA"
  set RHEAD = 0
else
  echo "Reported Rfactor     : $RFACT" | tee -a $LOG
  set RHEAD = 1
endif
if ($RFREE == 0.9990) then
  echo "Reported R-free      : none" | tee -a $LOG
  set RFHEAD = 0
  set RFREE = "NA"
else
  echo "Reported R-free      : $RFREE" | tee -a $LOG
  set RFHEAD = 1
endif
echo "Average B-factor     : $BAVER"    | tee -a $LOG
echo "Wilson B from SFcheck: $BWILS"    | tee -a $LOG
echo "Solvent percentage   : $SOLVD"    | tee -a $LOG
echo "Solvent model        : $SOLVENT"  | tee -a $LOG
echo "VdW probe            : $VDWPROBE" | tee -a $LOG
echo "Ion probe            : $IONPROBE" | tee -a $LOG
echo "Shrinkage            : $RSHRINK"  | tee -a $LOG
echo "Reflections (data)   : $REFCNT"   | tee -a $LOG
echo "Work set size        : $PWORKCNT" | tee -a $LOG
echo "Test set size (data) : $TSTCNT ($TSTPRC%)" | tee -a $LOG
echo "Test set size (used) : $PNTSTCNT " | tee -a $LOG
echo "Completeness (pdb)   : $COMPLETEH" | tee -a $LOG
echo "Completeness (used)  : $COMPLETED" | tee -a $LOG
echo "Wavelength           : $WAVELENGTH"| tee -a $LOG
if ($PTWINA != "") then
  echo "Twin fraction alpha  : $PTWINA"  | tee -a $LOG
endif
echo "Refinement tool      : $PROG"      | tee -a $LOG
echo "Deposition year      : $DYEAR"     | tee -a $LOG


############################################# Recalculation of R and R-free ##############################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG

#Stop if the completeness is too low (less than 20% or less than two thirds of reported completeness)
if (`perl -e 'if ($ENV{COMPLETED} < 20.0) {if (abs($ENV{COMPLETED} - $ENV{PCOMPLETEH}) < 2.5) {print "0"} else {print "1"}} elsif ($ENV{COMPLETED} < 0.67*$ENV{PCOMPLETEH}) {print "1"} else {print "0"}'` == 1) then
  echo "-This data is much too incomplete to use" | tee -a $LOG
  echo " o Cannot continue" | tee -a $LOG
  echo "COMMENT: Too much missing experimental data" >> $WHYNOT
  echo "PDB-REDO,$PDBID"                             >> $WHYNOT
  if ($SERVER == 1) then
    #Write out status files
    touch $STDIR/stoppingProcess.txt
    touch $STDIR/processStopped.txt
  endif
  cd $BASE
  exit(1)
endif

echo "****** R(-free) recalculation ******" | tee -a $LOG

#NCS settings
if ($DONCS == 1) then
  set NCSTYPE  = "ncsr local"
  set NCSALIGN = "ncsr align level 0.90 iterate N rmslevel 2.00"
  set NCSNEIGH = "ncsr neighbours exclude"
endif

#Setup Refmac library file for ligands
if ("$INREST" != "") then
  echo "-Importing geometric restraints" | tee -a $LOG
  cp $INREST $WORKDIR/${PDBID}_het.cif
  #Test the integrety of the restraint file
  if (`grep -a -c '_chem_comp_bond.value_dist' $WORKDIR/${PDBID}_het.cif` == 0 && `grep -a -c '_chem_link_bond.value_dist' $WORKDIR/${PDBID}_het.cif` == 0) then
    #There are no distance restraints. Report...
    echo " o The restraint file contained no valid restraints"   | tee -a $LOG
    echo " o Using standard and on-the-fly generated restraints" | tee -a $LOG
    echo "COMMENT: Corrupt user-provided restraints" >> $DEBUG
    echo "PDB-REDO,$PDBID"                           >> $DEBUG
    # ...and use the regular restraint generation
    set INREST = ""
    set LIBLIN = `echo LIBOUT $WORKDIR/${PDBID}_het.cif` #Output Refmac library file for new compounds or LINKs
  else
    #Force Refmac to use the uploaded restraint file
    set LIBLIN = `echo LIBIN $WORKDIR/${PDBID}_het.cif LIBOUT $WORKDIR/${PDBID}_het2.cif`
  endif
else
  set LIBLIN = `echo LIBOUT $WORKDIR/${PDBID}_het.cif` #Output Refmac library file for new compounds or LINKs
endif

#Import the additional restraints
if ("$INEXT" != "") then
  echo "-Importing additional restraints" | tee -a $LOG
  #Check whether the file is not a cif restraint file or a PDB file or a TLS definition
  if (`grep -a -c '_chem_comp_bond.value_dist' $INEXT` == 0 && `grep -a -c '_chem_link_bond.value_dist' $INEXT` == 0 && `grep -a -c '^[AH][TE][OT][MA]' $INEXT` == 0 && `grep -a -c '^RANGE' $INEXT` == 0 && `grep -c 'exte' $INEXT` != 0) then
    #Set up the external restraint file
    echo "external UndefinedAtoms ignore"  > $WORKDIR/external.rest
    echo "external weight scale $EXTSCAL" >> $WORKDIR/external.rest
    sed 's/  / /g' $INEXT >> $WORKDIR/external.rest
    #Set-up Refmac to use the uploaded additional restraints
    set RESTCMD = "@$WORKDIR/external.rest"
  else
    #This is not the right type of restraint file. Report...
    echo " o The additional restraint file is not in REFMAC format; it will be ignored" | tee -a $LOG
    echo "COMMENT: Corrupt external restraints" >> $DEBUG
    echo "PDB-REDO,$PDBID"                      >> $DEBUG
  endif
endif


#Set up scaling function
set SCALING = "lssc function a sigma $EXPSIG"

#Restarting oint with previous PDB-REDO entry
previousredo:

#TLS settings
junkedtls:
set ORITLS = 0
set ISTLS  = 'notls'
set TLSLIN =     #Do not use static TLS tensors unless...
#.... a TLS exists...
if (-e $WORKDIR/$PDBID.tls && $DOTLS == 1) then
  #...and has complete TLS tensors.
  if (`grep -c ORIGIN $WORKDIR/$PDBID.tls` != 0) then
    set ORITLS = 1
    set ISTLS  = 'tls'
    set TLSLIN = `echo TLSIN $WORKDIR/$PDBID.tls` #Use static TLS-groups for better R/R-free calculation
    echo -n "-Recalculating R-factors with original TLS tensors    " | tee -a $LOG
  else
    echo -n "-Running Refmac for 0 cycles "    | tee -a $LOG
  endif
else
  echo -n "-Running Refmac for 0 cycles " | tee -a $LOG
endif

#Run refmac for 0 cycles to check R-factors (with TLS). PROGRAM: REFMAC
refmac5 \
XYZIN  $WORKDIR/pdb${PDBID}.ent \
XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
HKLIN  $WORKDIR/$PDBID.mtz \
HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
$LIBLIN \
$TLSLIN \
$SCATLIN \
<<eof >& $WORKDIR/${PDBID}_0cyc$ISTLS.log
  $SCATTERCMD
  make check NONE
  make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
    ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
  refi type REST resi MLKF meth CGMAT bref MIXE
  $REFIRES
  ncyc 0
  tlsd waters exclude
  scal type $SOLVENT $SCALING
  solvent YES
  $MASKPAR
  $LOWMEM
  weight $WGTSIG MATRIX 0.5
  monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
    chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
  $NCSTYPE
  $NCSALIGN
  $NCSNEIGH
  $NCSSTRICT
  labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
  $ANOMCMD
  pdbout copy remarks 280 350
  NOHARVEST
  END
eof

if ($status || `grep -c 'Error: Fatal error. Cannot continue' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
  echo " " | tee -a $LOG
  #Check whether there are too many TLS goups
  if (`grep -c 'Too many tls groups. Maximum allowed is' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
    #Give error message
    echo " o Input model had more TLS groups than Refmac can handle. They will be ignored." | tee -a $LOG
    echo "COMMENT: too many TLS groups" >> $DEBUG
    echo "PDB-REDO,$PDBID"              >> $DEBUG  
    
    #Recover by ignoring te original TLS model
    mv $WORKDIR/$PDBID.tls $WORKDIR/$PDBID.notls
    goto junkedtls
  endif
  #Check to see if there is a problem with alternate residues
  if (`grep -c 'different residues have the same number' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
    if ($LOCAL == 1) then
      #Give the long error message
      echo " " | tee -a $LOG
      echo "FATAL ERROR!" | tee -a $LOG
      echo "------------" | tee -a $LOG
      echo "Refmac had problems using these residues with alternate identities:" | tee -a $LOG
      grep 'ERROR:' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 12- | tee -a $LOG
      echo "This problem may be solved by renumbering residues or by ensuring that" | tee -a $LOG
      echo "alternate atoms directly follow eachother in the PDB file." | tee -a $LOG
      echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details."  | tee -a $LOG
    else
      #Give the simple error message
      echo " o Cannot use structure with alternate residues" | tee -a $LOG
    endif
    #Write out WHY_NOT mesage
    echo "COMMENT: refmac: error with alternate residues" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
  #Check to see if there are dictionary problems
  if (! -e $WORKDIR/${PDBID}_het.cif || "$INREST" != "") then
    #Check for atoms not in the dictionary
    if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then

      #Give PDB_REDO mode-specific error messages
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Some residues have (atoms with) naming conflicts." | tee -a $LOG
        echo "Please, use standard atom and residue names in your input PDB file or upload a custom restraint file." | tee -a $LOG
        echo " " | tee -a $LOG
        echo "Residue  PDB standard description" | tee -a $LOG
        foreach HETID (`grep 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 22-37 | sort -u`)
          set HETID1 = `echo $HETID |cut -c 1-1`
          #Is the entry in LigandExpo
          wget --spider -q http://ligand-expo.rcsb.org/reports/$HETID1/$HETID
          if ($status) then
            echo "$HETID      New compound: make sure all $HETID residues are consistent within the input PDB" | tee -a $LOG
          else
            echo "$HETID      http://ligand-expo.rcsb.org/reports/$HETID1/$HETID" | tee -a $LOG
          endif
        end

        #Give server-specific extra information
        if ($SERVER == 1) then
          echo "Details from Refmac output:" | tee -a $LOG
          grep ' ERROR :' $WORKDIR/${PDBID}_0cyc$ISTLS.log | tee -a $LOG
        else
          echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details." | tee -a $LOG
        endif
      else
        #Give the short error message
        echo " o Residue or atom naming conflict. Cannot continue." | tee -a $LOG
      endif

      #Give WHY_NOT comment
      echo "COMMENT: refmac: residue or atom naming conflict" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                  >> $WHYNOT

    #Check for problems making a compound description
    else if (`grep -a -c 'is not completely connected' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then

      #Give PDB_RED mode-specific error messages
      if ($LOCAL == 1) then
        #Give the long error message
        echo " " | tee -a $LOG
        echo "FATAL ERROR!" | tee -a $LOG
        echo "------------" | tee -a $LOG
        echo "Cannot create restraint file for residue:" | tee -a $LOG
        grep 'program will create complete description for' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -d ':' -f 2 | sort -u | tee -a $LOG
        echo "Please, supply a restraint file using '--restin=Your_restraints.cif'." | tee -a $LOG
      else
        #Give the short error message
        echo " o Cannot create restraint file. Cannot continue." | tee -a $LOG
      endif

      #Write WHY_NOT comment
      echo "COMMENT: refmac: cannot create restraint file" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                               >> $WHYNOT

    #Check for NCS alignment problems
    else if (`grep -a -c 'ncs_ncs_generate.f90' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0 && ! -e $WORKDIR/renumber.json) then
        
      #Renumber the PDB file and start again
      cp $WORKDIR/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.bak
      $TOOLS/rnbterror.py -v -j $WORKDIR/renumber.json $WORKDIR/pdb${PDBID}.bak $WORKDIR/pdb${PDBID}.ent > $WORKDIR/renumber.log
          
      #Write DEBUG message
      echo " o Problem with NCS alignment. Renumbering terminal residues and restarting." | tee -a $LOG
      echo "COMMENT: residues renumbered" >> $DEBUG
      echo "PDB-REDO,$PDBID"              >> $DEBUG  
     
     #Start again
      goto renumbered

    #All other problems
    else
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
    endif

    #Write out status files
    if ($SERVER == 1) then
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif


#Check to see whether a new ligand was encountered (only if no extra restraints were provided).
if ("$INREST" == "") then
  if(-e $WORKDIR/${PDBID}_het.cif) then
    echo " " | tee -a $LOG
    echo -n " o New ligand encountered, retrying " | tee -a $LOG

    #Backup old files
    cp $WORKDIR/${PDBID}_0cyc$ISTLS.log $WORKDIR/${PDBID}_0cycv1.log
    if(-e $WORKDIR/${PDBID}_0cyc$ISTLS.pdb) then
      cp $WORKDIR/${PDBID}_0cyc$ISTLS.pdb $WORKDIR/${PDBID}_0cycv1.pdb
    endif

    #Start using the new ligand library
    set LIBLIN = `echo LIBIN $WORKDIR/${PDBID}_het.cif LIBOUT $WORKDIR/${PDBID}_het2.cif`

    #Rerun refmac for 0 cycles
    refmac5 \
    XYZIN  $WORKDIR/pdb${PDBID}.ent \
    XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
    $LIBLIN \
    $TLSLIN \
    $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_0cyc$ISTLS.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
	ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 280 350
      NOHARVEST
      END
eof
    if($status) then
      #Try to give a specific error message
      #Check to see if there is a problem with alternate residues
      if (`grep -c 'different residues have the same number' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
        if ($LOCAL == 1) then
          #Give the long error message
          echo " " | tee -a $LOG
          echo "FATAL ERROR!" | tee -a $LOG
          echo "------------" | tee -a $LOG
          echo "Refmac had problems using these residues with alternate identities:" | tee -a $LOG
          grep 'ERROR:' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 12- | tee -a $LOG
          echo "This problem may be solved by renumbering residues or by ensuring that" | tee -a $LOG
          echo "alternate atoms directly follow eachother in the PDB file." | tee -a $LOG
          echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details."  | tee -a $LOG
        else
          #Give the simple error message
          echo " " | tee -a $LOG
          echo " o Cannot use structure with alternate residues" | tee -a $LOG
        endif
        #Write out WHY_NOT mesage
        echo "COMMENT: refmac: error with alternate residues" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      endif

      #Check to see if there are dictionary problems
      if (! -e $WORKDIR/${PDBID}_het2.cif) then
        #Check for atoms not described in the dictionary
        if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
          if ($LOCAL == 1) then
            #Give the long error message
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Some residues have (atoms with) naming conflicts." | tee -a $LOG
            echo "Please, use standard atom and residue names in your input PDB file or upload a custom restraint file." | tee -a $LOG
            echo " " | tee -a $LOG
            echo "Residue  PDB standard description" | tee -a $LOG
            echo "---------------------------------" | tee -a $LOG
            foreach HETID (`grep 'is absent in the library' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -c 22-37 | sort -u`)
              set HETID1 = `echo $HETID |cut -c 1-1`
              #Is the entry in LigandExpo
              wget --spider -q http://ligand-expo.rcsb.org/reports/$HETID1/$HETID
              if ($status) then
                echo "$HETID      New compound: make sure all $HETID residues are consistent within the input PDB" | tee -a $LOG
              else
                echo "$HETID      http://ligand-expo.rcsb.org/reports/$HETID1/$HETID" | tee -a $LOG
              endif
            end
            echo " " | tee -a $LOG

            #Give server-specific extra information
            if ($SERVER == 1) then
              echo "Details from Refmac output:" | tee -a $LOG
              grep ' ERROR :' $WORKDIR/${PDBID}_0cyc$ISTLS.log | tee -a $LOG
            else
              echo "See the file $WORKDIR/${PDBID}_0cyc$ISTLS.log for details." | tee -a $LOG
            endif
          else
            #Give the short error message
            echo " " | tee -a $LOG
            echo " o Residue or atom naming conflict. Cannot continue." | tee -a $LOG
          endif
          echo "COMMENT: refmac: residue or atom naming conflict" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                  >> $WHYNOT
        #Check for residues for which a restraint file cannot be generated
        else if (`grep -a -c 'is not completely connected' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0) then
          if ($LOCAL == 1) then
            #Give the long error message
            echo " " | tee -a $LOG
            echo "FATAL ERROR!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Cannot create restraint file for residue:" | tee -a $LOG
            grep 'program will create complete description for' $WORKDIR/${PDBID}_0cyc$ISTLS.log | cut -d ':' -f 2 | sort -u
            echo "Please, supply a restraint file using '--restin=Your_restraints.cif'." | tee -a $LOG
          else
            #Give the short error message
            echo " " | tee -a $LOG
            echo " o Cannot create restraint file. Cannot continue." | tee -a $LOG
          endif
          echo "COMMENT: refmac: cannot create restraint file" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                               >> $WHYNOT
        #Check for problems with local NCS restraints
        else if (`grep -a -c 'ncs_ncs_generate.f90' $WORKDIR/${PDBID}_0cyc$ISTLS.log` != 0 && ! -e $WORKDIR/renumber.json) then
          #Renumber the PDB file and start again
          cp $WORKDIR/pdb${PDBID}.ent $WORKDIR/pdb${PDBID}.bak
          $TOOLS/rnbterror.py -j $WORKDIR/renumber.json $WORKDIR/pdb${PDBID}.bak $WORKDIR/pdb${PDBID}.ent > $WORKDIR/renumber.log
          
          #Write DEBUG message
          echo " " | tee -a $LOG
          echo " o Problem with NCS alignment. Renumbering terminal residues and restarting." | tee -a $LOG
          echo "COMMENT: residues renumbered" >> $DEBUG
          echo "PDB-REDO,$PDBID"              >> $DEBUG  
          
          #Start again
          goto renumbered
         
        else
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
        endif
      endif
      #Stop the run
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  endif
endif

#Check the value for LOGSTEP and RBLS
if (`grep -a -c 'Rms ChirVolume' $WORKDIR/${PDBID}_0cyc$ISTLS.log` == 0) then
  @ LOGSTEP = ($LOGSTEP - 1)
  @ RBLS = ($LOGSTEP - 2)
else
  @ RBLS = ($LOGSTEP - 3)
endif

#Get calculated R-factor
setenv PRCAL1  `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $2}'`
set    TRFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $3}'`
echo "(R-free = $TRFREE)" | tee -a $LOG

#Check R-factor without TLS tensors
if ($ORITLS == 1) then
  echo -n "-Recalculating R-factors without original TLS tensors " | tee -a $LOG
  set ISTLS = 'notls'

  #Rerun refmac for 0 cycles, again
  refmac5 \
  XYZIN  $WORKDIR/pdb${PDBID}.ent \
  XYZOUT $WORKDIR/${PDBID}_0cyc$ISTLS.pdb \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_0cyc$ISTLS.mtz \
  $LIBLIN \
  $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc$ISTLS.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type REST resi MLKF meth CGMAT bref MIXE
    $REFIRES
    ncyc 0
    tlsd waters exclude
    scal type $SOLVENT $SCALING
    solvent YES
    $MASKPAR
    $LOWMEM
    weight $WGTSIG MATRIX 0.5
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $NCSTYPE
    $NCSALIGN
    $NCSNEIGH
    $NCSSTRICT
    labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
    $ANOMCMD
    pdbout copy remarks 280 350
    NOHARVEST
    END
eof
  if ($status) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Get calculated R-factor
  setenv PRCAL2  `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $2}'`
  set    TRFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$ISTLS.log | head -n 1 | awk '{print $3}'`
  echo "(R-free = $TRFREE)" | tee -a $LOG

  #Which is better (i.e. gives lower R-factor), with or without TLS?
  if (`perl -e 'if ($ENV{PRCAL1} - $ENV{PRCAL2} < 0) {print "with"} else {print "without"};'` == 'with') then
    echo " o Recalculating the R-factor with TLS tensors gives the best result" | tee -a $LOG
    set ISTLS = 'tls'
  else
    echo " o Recalculating the R-factor without TLS tensors gives the best result" | tee -a $LOG
    set TLSLIN =   #Do not use static TLS tensors
  endif
endif

#Clean up a bit
mv $WORKDIR/${PDBID}_0cyc$ISTLS.log $WORKDIR/${PDBID}_0cyc.log
mv $WORKDIR/${PDBID}_0cyc$ISTLS.pdb $WORKDIR/${PDBID}_0cyc.pdb
mv $WORKDIR/${PDBID}_0cyc$ISTLS.mtz $WORKDIR/${PDBID}_0cyc.mtz

#Evaluate log file (dirty)
set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`

#Check whether de-twinning is needed to repoduce R-factors. Only when NOT using legacy mode
if ($LEGACY == 0 && $DOTWIN == 1) then
  #Did SFCHECK report a possible twin
  set TWIN = `perl -e 'if ($ENV{PTWINA} < 0.05) {print ""} else {print "test"};'`

  #Check whether the R-factors can be reproduced
  if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0 && $TWIN == 'test') then #R-factors do not fit, try de-twinning.
    echo " o Problem reproducing the R-factors" | tee -a $LOG
    echo -n "-Trying de-twinning " | tee -a $LOG

    #Set up detwinning
    set TWIN = 'twin'

    #Create a backup
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv3.log

    #Run Refmac
    refmac5 \
    XYZIN  $WORKDIR/pdb${PDBID}.ent \
    XYZOUT $WORKDIR/${PDBID}_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $TLSLIN \
    $LIBLIN \
    $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 280 350
      NOHARVEST
      END
eof
    if ($status) then
      echo " " | tee -a $LOG
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Switch off detwinning if Refmac finds just one twin domain
    if (`grep -a 'The number of twin domains' $WORKDIR/${PDBID}_0cyc.log | tail -n 1 | awk '{print $7}'` == 1) then
      set TWIN =  #No detwinning
      echo " " | tee -a $LOG
      echo " o Refmac detected no twinning" | tee -a $LOG
    else
      #Report R-free
      set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`
      echo "(R-free = $RFCAL)" | tee -a $LOG

      #Did Refmac find higher symmetry
      set REFHSYMM = `grep -c 'twin or higher symmetry' $WORKDIR/${PDBID}_0cyc.log`

      #Evaluate R-factors
      set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
      if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 1) then
        #Detwinning was needed to reproduce R-factors
        #Give warnings
        if ($PHASTWIN == 0 && $REFHSYMM > 0) then
          if ($LOCAL == 1) then
            #Give big warning for server and local mode
            echo " " | tee -a $LOG
            echo "WARNING!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Data seems to have treated as twinnned before, but Phaser did not find significant twinning." | tee -a $LOG
            echo "Additionally, Refmac finds a potential space group problem." | tee -a $LOG
            echo "The data will not be treated as twinned in PDB-REDO." | tee -a $LOG
            echo "Please, check the space group carefully!" | tee -a $LOG
            echo " " | tee -a $LOG

            #Do not use twinning and copy back the old log file
            cp $WORKDIR/${PDBID}_0cycv3.log $WORKDIR/${PDBID}_0cyc.log
            set TWIN =  #No detwinning

          else
            #Give small warning for databank mode
            echo " o Data seems to have treated as twinnned before."  | tee -a $LOG
            echo "   * Phaser did not find significant twinning."     | tee -a $LOG
            echo "   * Refmac finds a potential space group problem." | tee -a $LOG
            echo "   * Reluctantly treating the data as twinned."     | tee -a $LOG

            set FALSETWIN = 1

            #Write DEBUG comment
            echo "COMMENT: false twinning problem" >> $DEBUG
            echo "PDB-REDO,$PDBID"                 >> $DEBUG

            #The data are treated as twinned. Do not calculate anomalous maps
            set ANOMCOEF =
            set ANOMCMD  =
          endif
        else if ($PHASTWIN == 0) then
          if ($LOCAL == 1) then
            #Give big warning for server and local mode
            echo " " | tee -a $LOG
            echo "WARNING!" | tee -a $LOG
            echo "------------" | tee -a $LOG
            echo "Data seems to have treated as twinnned before and Refmac and sfcheck indicate that the data are twinned, but Phaser does not. " | tee -a $LOG
            echo "The data will be treated as twinned in PDB-REDO, reluctantly." | tee -a $LOG
            echo "Please, check your data carefully!" | tee -a $LOG
            echo " " | tee -a $LOG
          else
            #Give small warning for databank mode
            echo " o Data seems to have treated as twinnned before." | tee -a $LOG
            echo "   * Phaser did not find significant twinning."    | tee -a $LOG
            echo "   * Reluctantly treating the data as twinned."    | tee -a $LOG

            set FALSETWIN = 1
            #Write DEBUG comment
            echo "COMMENT: ambiguous twinning problem" >> $DEBUG
            echo "PDB-REDO,$PDBID"                     >> $DEBUG
          endif

          #The data are treated as twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        else
          #The data are twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        endif
      else
        #Detwinning did not solve the problem. Only use twinning when it is detected
        if ($PHASTWIN == 1 && $REFHSYMM == 0) then
          #The data are twinned. Do not calculate anomalous maps
          set ANOMCOEF =
          set ANOMCMD  =
        else
          #Do not use twinning and copy back the old log file
          cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv4.log
          cp $WORKDIR/${PDBID}_0cycv3.log $WORKDIR/${PDBID}_0cyc.log
          set TWIN =  #No detwinning
        endif

      endif
    endif
  endif

  #Evaluate log file (dirty)
  set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
  set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`

endif


#Check the fit with the data. If it is poor, try rigid-body refinement. This is always done when in legacy mode
if ($DORB == 1 && (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0 || $LEGACY == 1)) then #R-factors do not fit, try rigid-body refinement.

  #Create a backup
  cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv5.log

  #Do $RBCYCLE cycles of rigid body refinement
  if ($LEGACY == 1) then
    echo "-Legacy mode" | tee -a $LOG
    echo -n " o Doing rigid-body refinement " | tee -a $LOG
    set RBCYCLE = 15
  else
    echo " o Problem reproducing the R-factors"  | tee -a $LOG
    echo -n "-Trying rigid-body refinement " | tee -a $LOG
  endif

  #Set up NCS
  if ($STRICTNCS == 1) then
    set RBNCS = ncsconstraints
  endif

  #Run Refmac
  refmac5 \
  XYZIN $WORKDIR/pdb$PDBID.ent \
  XYZOUT $WORKDIR/${PDBID}_refmacrb.pdb \
  HKLIN $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_refmacrb.mtz \
  $SCATLIN \
<<eof > $WORKDIR/${PDBID}_rb.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type RIGID resi MLKF meth CGMAT
    $REFIRES
    mode rigid
    rigid ncycle $RBCYCLE
    scal type $SOLVENT $SCALING
    scale mlscale nrfr 5
    solvent YES
    $MASKPAR
    $LOWMEM
    weight $WGTSIG MATRIX 0.5
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $RBNCS
    $TWIN
    labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
    $ANOMCMD
    pdbout copy remarks 280 350
    NOHARVEST
    END
eof
  if ($status) then
    echo " " | tee -a $LOG
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in rigid-body refinement" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Transplant the LINKs if needed
  if (`grep -c ^LINK $WORKDIR/pdb$PDBID.ent` != 0 && `grep -c ^LINK $WORKDIR/${PDBID}_refmacrb.pdb` == 0) then
    cp  $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_refmacrb.bak
    grep -B 100000 ^CRYST1 $WORKDIR/${PDBID}_refmacrb.bak | grep -v ^CRYST1 > $WORKDIR/${PDBID}_refmacrb.pdb
    grep ^LINK $WORKDIR/pdb$PDBID.ent >> $WORKDIR/${PDBID}_refmacrb.pdb
    grep -A 250000 ^CRYST1 $WORKDIR/${PDBID}_refmacrb.bak >> $WORKDIR/${PDBID}_refmacrb.pdb
  endif


  #Now do another restrained refinement run to include the TLS contribution (if needed)
  if ($ORITLS == 1) then
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refmacrb.pdb \
    XYZOUT $WORKDIR/${PDBID}_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $TLSLIN \
    $LIBLIN \
    $SCATLIN \
<<eof > $WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref MIXE
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 280 350
      NOHARVEST
      END
eof
    if ($status) then
      echo " " | tee -a $LOG
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Evaluate log file
    set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
    set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`
  else
    #Evaluate log file from the rigid body refinement
    set RCAL  = `tail -n $RBLS $WORKDIR/${PDBID}_rb.log | head -n 1 | awk '{print $2}'`
    set RFCAL = `tail -n $RBLS $WORKDIR/${PDBID}_rb.log | head -n 1 | awk '{print $3}'`

    #Use the MTZ file from the rigid-body refinement
    cp $WORKDIR/${PDBID}_refmacrb.mtz $WORKDIR/${PDBID}_0cyc.mtz
  endif

  #Fill the line with the R-free value
  echo "(R-free = $RFCAL)" | tee -a $LOG

  #Copy files so that the rigid-body refined model is used
  cp $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_0cyc.pdb
  
else  
  cp $WORKDIR/pdb$PDBID.ent $WORKDIR/${PDBID}_refmacrb.pdb
endif

#Check the fit with the data. If it is poor, try short TLS-refinement
if (`$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL` == 0) then
  #Only do something if original TLS groups were given
  if ($ORITLS == 1) then
    echo " o Problem reproducing the R-factors" | tee -a $LOG
    echo -n "-Trying TLS tensor re-evaluation " | tee -a $LOG

    #Create a backup
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cycv5.log
    cp $WORKDIR/${PDBID}_0cyc.mtz $WORKDIR/${PDBID}_0cycv5.mtz

    #Use TLS group definitions only
    grep -E 'TLS|RANGE|^ ' $WORKDIR/${PDBID}.tls > $WORKDIR/${PDBID}_0cycin.tls

    #Run Refmac with 5 TLS cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refmacrb.pdb \
    XYZOUT $WORKDIR/${PDBID}_TLS0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_0cyc.mtz \
    $LIBLIN \
    TLSIN $WORKDIR/${PDBID}_0cycin.tls TLSOUT $WORKDIR/${PDBID}_0cyc.tls \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_0cyc.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi tlsc 5
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc 0
      tlsd waters exclude
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG MATRIX 0.5
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $TWIN
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES  $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 280 350
      NOHARVEST
      kill $TOOLS/pdb_redo.refmac
      END
eof
    if ($status) then
      if (`grep -a -c 'Program terminated by user' $WORKDIR/${PDBID}_0cyc.log` != 0) then
        #Problems with the TLS group definition.
        mv $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.notls
        #Cannot use this run. Resore the backup.
        cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      else
        echo " " | tee -a $LOG
        echo " o Problem with refmac. Cannot continue." | tee -a $LOG
        echo "COMMENT: refmac: error in initial R-free calculation" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      endif
    else if (`grep -c -E '\*{6}' $WORKDIR/${PDBID}_TLS0cyc.pdb` != 0) then
      #TLS refinement caused rediculously high B-factors. Don't use this TLS group selection.
      mv $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.notls
      #Cannot use this run. Resore the backup.
      cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      #Report
      echo " " | tee -a $LOG
      echo " o TLS refinement was unstable. Original TLS group selection will not be used." | tee -a $LOG
      echo "COMMENT: Problem with TLS tensor re-evaluation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                >> $DEBUG
    else
      #Strip ANISOU records and clean up
      grep -a -v -E '^ANISOU' $WORKDIR/${PDBID}_TLS0cyc.pdb > $WORKDIR/${PDBID}_0cyc.pdb
      rm $WORKDIR/${PDBID}_TLS0cyc.pdb
      rm $WORKDIR/${PDBID}_0cycin.tls
      rm $WORKDIR/${PDBID}_0cyc.tls

      #Did the R-factor improve over the rigid-body results?
      setenv PRCAL1  `echo $RCAL`
      setenv PRCAL2  `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`

      #Discard the TLS results if they are poor.
      if (`perl -e 'if ($ENV{PRCAL1} - $ENV{PRCAL2} < 0) {print "rb"} else {print "tls"};'` == 'rb') then
        cp $WORKDIR/${PDBID}_refmacrb.pdb $WORKDIR/${PDBID}_0cyc.pdb
        cp $WORKDIR/${PDBID}_0cycv5.mtz $WORKDIR/${PDBID}_0cyc.mtz
        cp $WORKDIR/${PDBID}_0cycv5.log $WORKDIR/${PDBID}_0cyc.log
      endif

      #Evaluate log file (dirty)
      set RCAL  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $2}'`
      set RFCAL = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | awk '{print $3}'`

      echo "(R-free = $RFCAL)" | tee -a $LOG
    endif
  endif
endif

#Get the geometry
set RMSZB = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $8}'`
set RMSZA = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc.log | head -n 1 | sed 's/\*\*\*\*\*\*/huge/g' | awk '{print $10}'`
setenv PRMSZB $RMSZB
setenv PRMSZA $RMSZA

#Check if R-factor can (very roughly) be reproduced. If not, abort. Do not do this in legacy mode.
if ($LEGACY == 1) then
  #Write a warning if R-free is really high
  if (`echo $RFCAL | awk '{if ($1 > 0.500) {print 1} else {print 0}}'` == 1) then
    echo "COMMENT: Suspiciously high calculated R-free" >> $DEBUG
    echo "PDB-REDO,$PDBID"                              >> $DEBUG
  endif
else
  set FITR = `$TOOLS/fitr $RFACT $RFREE $RCAL $RFCAL 0.10`
  if ($FITR == 0) then #R-factors do not fit after several tries.
    echo " o Cannot reproduce Rfactor within 0.10 tolerance:" | tee -a $LOG
    echo "   * Reported Rfactor : $RFACT" | tee -a $LOG
    echo "   * Calculated R     : $RCAL"  | tee -a $LOG

    if ($INTENS != "-i" && `grep -c ' F ' $WORKDIR/${PDBID}c2c.log` > 0 &&  `grep -c ' I ' $WORKDIR/${PDBID}c2c.log` > 0) then

      #Go back and use intensities
      echo "   * Retrying with reflection intensities intead of amplitudes"  | tee -a $LOG
      echo " " | tee -a $LOG

      #Create debug entry
      echo "COMMENT: Using intensities instead of amplitudes" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                  >> $DEBUG

      #Set the flag for cif2cif and go back there
      set INTENS = "-i"
      set SIGMA = ''
      goto beintens
    else if ($LOCAL == 0 && $GOTOLD == 0) then
      #See if a previous PDB-REDO entry exists and use its 0-cycle file
      $WEBGET http://www.cmbi.ru.nl/pdb_redo/$D2/$PDBID/${PDBID}_0cyc.pdb.gz
      if($status) then
        echo "   * Re-refinement aborted" | tee -a $LOG
      else
        zcat $WORKDIR/${PDBID}_0cyc.pdb.gz > $WORKDIR/pdb$PDBID.ent
        echo "   * Trying previous PDB-REDO entry as starting point" | tee -a $LOG
        set GOTOLD = 1
        goto previousredo
      endif
    else
      echo "   * Re-refinement aborted" | tee -a $LOG
    endif

    #Create whynot entry
    echo "COMMENT: Cannot reproduce Rfactor within 0.10 tolerance" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
endif

#Check whether Refmac made anomalous maps
if (`grep -c 'No anomalous scatterers found' $WORKDIR/${PDBID}_0cyc.log` > 0) then
  #Refmac will not make anomalous maps
  set ANOMCOEF =
  set ANOMCMD  =
endif


####################################### Decide on using detwinning in refinement #########################################

#If detwinning is active, test to see if Refmac detects a twin
if ($TWIN == 'test' && $PHASTWIN == 1) then
  echo "-Evaluating twinning" | tee -a $LOG

  #Set up detwinning
  set TWIN = 'twin'

  refmac5 \
  XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
  XYZOUT $WORKDIR/${PDBID}_twin.pdb \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/${PDBID}_twin.mtz \
  $LIBLIN \
  $SCATLIN \
<<eof >$WORKDIR/${PDBID}_twin.log
    $SCATTERCMD
    make check NONE
    make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
      ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
    refi type REST resi MLKF meth CGMAT bref MIXE
    $REFIRES
    ncyc 0
    scal type $SOLVENT $SCALING
    solvent YES
    $MASKPAR
    $LOWMEM
    $NCSSTRICT
    weight $WGTSIG AUTO 2.50
    monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
      chiral  10.0 bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
    $TWIN
    labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES
    pdbout copy remarks 280 350
    END
eof
  if ($status) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in twin evaluation" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                           >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif

  #Switch off detwinning if Refmac finds just one twin domain
  if (`grep -a 'The number of twin domains' $WORKDIR/${PDBID}_twin.log | tail -n 1 | awk '{print $7}'` == 1) then
    set TWIN =  #No detwinning
    echo " o Refmac detected no twinning" | tee -a $LOG
  else
    #The data are twinned. Do not calculate anomalous maps.
    set ANOMCOEF =
    set ANOMCMD  =
  endif

  #Clean up
  rm $WORKDIR/${PDBID}_twin.mtz
  rm $WORKDIR/${PDBID}_twin.pdb
else if ($TWIN == 'test') then
  #Do not try detwinning because PHASER detected no twinning
  set TWIN =  #No detwinning
endif

#Set the final twin status
if ($TWIN == "twin") then
  set ISTWIN = 1
endif

############################################## Check R-factors and geometry ##############################################

#Warn for high structure RMSZ scores
if ($RMSZB == 'huge' || $RMSZA == 'huge' || `perl -e 'if ($ENV{PRMSZB} > 10.000) {print "1"} elsif ($ENV{PRMSZA} > 10.000) {print "1"} else {print "0"};'` == 1 ) then
  if ($LOCAL == 1) then
    #The warning
    echo " " | tee -a $LOG
    echo "WARNING!" | tee -a $LOG
    echo "--------" | tee -a $LOG
    echo "Extremely large geometric outliers!" | tee -a $LOG
    echo "Bond length RMSZ: $RMSZB" | tee -a $LOG
    echo "Bond angle RMSZ : $RMSZA" | tee -a $LOG
    echo "This is likely a problem with the restraint generation." | tee -a $LOG
    if ("$INREST" != "") then
      echo "Make sure that all non-standard compounds and LINKs are properly described in your restraint file."   | tee -a $LOG
    else
      echo "Consider running PDB-REDO with a restraint file that describes all non-standard compounds and LINKs." | tee -a $LOG
    endif
    echo " " | tee -a $LOG
  else
    echo "-Suspiciously high Refmac bond length or bond angle RMSZ detected" | tee -a $LOG
    echo " o Bond length RMSZ: $RMSZB"                               | tee -a $LOG
    echo " o Bond angle RMSZ : $RMSZA"                               | tee -a $LOG
    echo " o This is likely a problem with the restraint generation" | tee -a $LOG
  endif

  #Write debug message
  echo "COMMENT: suspiciously high rmsZ" >> $DEBUG
  echo "PDB-REDO,$PDBID"                 >> $DEBUG
endif


#Set value for bits of Perl code in the script
setenv PRCAL  $RCAL
setenv PRFCAL $RFCAL
setenv PPPATM  `echo "4"` #Assumes isotropic B-factors for now

#Calculate sigma(R-free)
set SIGRFCAL = `perl -e 'printf ("%.4f\n", $ENV{PRFCAL}/sqrt($ENV{PNTSTCNT}));'`
setenv PSIGRFCAL $SIGRFCAL

#Calculate unbiased R-free/R ratio (thanks to Ian Tickle)
set RFRRAT = `perl -e '$x = $ENV{PATMCNT}/$ENV{PWORKCNT}; $z = 1-$ENV{PPPATM}*$x/(1+2.5*$x); $z = 0 if $z<0; $a = $ENV{PPPATM}-2.5*(1-$z**5); $x = $a*$x; $rfree_rat = 1-$x; if($rfree_rat <= 0) {$rfree_rat = 1.010} else {$rfree_rat = sqrt((1+$x)/$rfree_rat)} $rfree_rat = 1.200000 if ($ENV{URESO}>2.65 and $rfree_rat>1.200000); $rfree_rat = 1.200000 if ($ENV{URESO}>3.0 and $rfree_rat<1.011); $rfree_rat = 1.450000 if $rfree_rat>1.454; printf ("%.6f\n", $rfree_rat);'`
setenv PRFRRAT $RFRRAT

#Calculate unbiased R-free
set RFCALUNB = `perl -e 'printf ("%.4f\n", $ENV{PRCAL}*$ENV{PRFRRAT});'`
setenv PRFCALUNB $RFCALUNB

#R-free Z-score (comparison of R-free with its expected value)
set RFCALZ = `perl -e 'printf ("%.2f\n", ($ENV{PRFCALUNB}-$ENV{PRFCAL})/$ENV{PSIGRFCAL});'`
setenv PRFCALZ $RFCALZ

#Print values
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** R-factor and R-free details ******" | tee -a $LOG
echo "Calculated R     : $RCAL"     | tee -a $LOG
echo "Calculated R-free: $RFCAL"    | tee -a $LOG
echo "'Unbiased' R-free: $RFCALUNB" | tee -a $LOG
echo "sigma(R-free)    : $SIGRFCAL" | tee -a $LOG
echo "R-free Z-score   : $RFCALZ"   | tee -a $LOG


############################################ Fix atom chiralities if needed ##############################################

#Count the chirality problems
set CHIRERR = `grep -a -A 100 "Chiral volume deviations from the" $WORKDIR/${PDBID}_0cyc.log | grep -E "^.{19}mod" | wc -l`

#Are fixes needed?
if ($CHIRERR != 0) then

  #Report
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Chirality validation ******" | tee -a $LOG
  echo "-Found $CHIRERR chirality problems" | tee -a $LOG
  echo " o Running chiron"   | tee -a $LOG

  #Do the fixes and report
  mv $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_0cyc.old

  #PROGRAM: chiron
  $TOOLS/chiron -v $TOOLS/pdb_redo.dat $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_0cyc.old $WORKDIR/${PDBID}_0cyc.pdb > $WORKDIR/chiron.log
  set CHIFIX = `grep "ters fixed" $WORKDIR/chiron.log | awk '{print $5}'`
  echo "   * `grep 'tested' $WORKDIR/chiron.log`"         | tee -a $LOG
  echo "   * `grep 'ters fixed' $WORKDIR/chiron.log`"     | tee -a $LOG
  echo "   * `grep 'not fixed' $WORKDIR/chiron.log`"      | tee -a $LOG
  echo "   * `grep 'Unknown chiral' $WORKDIR/chiron.log`" | tee -a $LOG

  #Make a debug record for unfixable problems
  if (`grep -c Unknown $WORKDIR/chiron.log` != 0 && `grep Unknown $WORKDIR/chiron.log | awk '{print $5}'` != 0) then
    if ($LOCAL != 1) then
      echo "$PDBID :"                                           >> $CHIRALS
      grep 'Unknown residue' $WORKDIR/chiron.log | cut -c 19-61 >> $CHIRALS
    else
      echo "COMMENT: Unknown chirality errors" >> $DEBUG
      echo "PDB-REDO,$PDBID"                   >> $DEBUG
    endif
  endif

  #Make a debug record if chiron failed
  if (-e $WORKDIR/${PDBID}_0cyc.pdb) then
    #Do nothing
  else
    mv $WORKDIR/${PDBID}_0cyc.old $WORKDIR/${PDBID}_0cyc.pdb
    echo "COMMENT: CHIRON: general error" >> $DEBUG
    echo "PDB-REDO,$PDBID"                >> $DEBUG
  endif
else
  set CHIFIX = 0
endif


###############################################  Set up occupancy fixing  ################################################
#Set counter
set OCCREF = 0

#Set up occupancy refinement if it is not surpressed
if ($DOOCC == 1) then
  #Create file with hetero compounds with more than two occupancies if needed
  if (`grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb |  grep -v MSE | cut -c 22-27,56-60 | sort -u | cut -c 1-6 | uniq -c | awk '{if ($1 > 2) {$1 = ""; print substr($0,2,6)}}' | wc -l` != 0) then
    #Grep the hetero compounds, but leave out MSE (seleno-methionine), then cut out the chains, residue numbers and
    #occupancies and sort. This will leave only the unique occupancies per residue. The chains and residue numbers are
    #cut out and for each residue the count is given. If this count is > 2, the chain and residue number is printed to
    # a file.
    grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb | grep -v MSE | cut -c 22-27,56-60 | sort -u | cut -c 1-6 | uniq -c | awk '{if ($1 > 2) {print substr($0,length($0)-5, length($0))}}' > $WORKDIR/occupancy.lst

  endif

  #Find hetero compounds with at least one atom with occupancy 0.01
  foreach RES ("`grep -E '^HETATM'+'.{10} ' $WORKDIR/${PDBID}_0cyc.pdb | grep -v MSE | cut -c 22-27,56-60 | grep ' 0.01' | cut -c 1-6 | sort -u`")
    echo "$RES" >> $WORKDIR/occupancy.lst
  end

  #Filter the residue list and write the Refmac command file
  if (-e $WORKDIR/occupancy.lst) then
    foreach RES ("`sort -u $WORKDIR/occupancy.lst`")
      @ OCCREF = ($OCCREF + 1)
      set CHNID = `echo $RES | cut -c 1-1`
      set RESN  = `echo $RES | cut -c 2-6`
      #Skip cases where the residue has an insertion code (i.e. $RESN has non-numeric characters)
      if (`echo $RESN | awk '{if ($0 ~ /^[0-9]+$/) {print "1"} else {print "0"}}'` == 1 ) then
        echo "occupancy group id $OCCREF chain $CHNID residue $RESN" >> $WORKDIR/occupancy_cmd.refmac
      else
        @ OCCREF = ($OCCREF - 1)
      endif
    end
  endif

  #Is occupancy refinement needed?
  if (-e $WORKDIR/occupancy_cmd.refmac) then

    #Append refinement command
    echo "occupancy refine" >> $WORKDIR/occupancy_cmd.refmac

    #Make Refmac use command file
    set OCCCMD = "@$WORKDIR/occupancy_cmd.refmac"
  endif
endif

############################################### Resolution-based settings ################################################
#Label to come back to after updating the resolution
resobasedsettings:

#Find the resolution type (6 categories, do not use reflections per atom for heavy strict NCS)
# Resolution categories: 0 = extremely low, 1 = very low, 2 = low, 3 = medium, 4 = high, 5 = atomic)
if ($STRICTNCS == 1 && $PMTRIX > 9) then
  set RESOTYPE = `perl -e 'if ($ENV{URESO} > 4.99) {print "0"} elsif ($ENV{URESO} > 3.49 and $ENV{URESO} < 5.00) {print "1"} elsif ($ENV{URESO} > 2.79 and $ENV{URESO} < 3.50) {print "2"} elsif ($ENV{URESO} < 1.21) {print "5"} elsif ($ENV{URESO} > 1.20 and $ENV{URESO} < 1.71) {print "4"} else {print "3"};'`
else
  set RESOTYPE = `perl -e 'if ($ENV{URESO} > 4.99) {print "0"} elsif ($ENV{PWORKCNT}/$ENV{PATMCNT} < 1.0) {print "0"} elsif ($ENV{URESO} > 3.49 and $ENV{URESO} < 5.00) {print "1"} elsif ($ENV{PWORKCNT}/$ENV{PATMCNT} < 2.5) {print "1"} elsif ($ENV{URESO} > 2.79 and $ENV{URESO} < 3.50) {print "2"} elsif ($ENV{URESO} < 1.21) {print "5"} elsif ($ENV{URESO} > 1.20 and $ENV{URESO} < 1.71) {print "4"} else {print "3"};'`
endif

#Modify the resolution category based on user input and experiment type
if ($ISED == 1) then
  @ RESOTYPE = ($RESOTYPE - 1)
endif
@ RESOTYPE = ($RESOTYPE + $RESCATSTEP)

#Set refinement parameters
if ($RESOTYPE < 1) then
  set WEIGHTS    = "1e-7 1e-6 1e-5 1e-4 5e-4 .001"
  set BWEIGHTS   = "2.50 2.00 1.50 1.20 1.00"  #Do not try looser B-factor restraints
  set JELLY      = "ridg dist sigm 0.02"       #Use jelly-body refinement
  @ NCYCLE = ($NCYCLE + 15)
  set TORSION    = "restr tors include group peptide"
  set ESTRICT    = "-e"                        #Picker is extra strict
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 1) then
  set WEIGHTS  = "5e-4 .001 .002 .005 0.01 0.03"
  set BWEIGHTS = "2.00 1.50 1.20 1.00 0.80 0.50"
  set JELLY    = "ridg dist sigm 0.05" #Use jelly-body refinement
  @ NCYCLE = ($NCYCLE + 10)
  set ESTRICT  = "-e"                  #Picker is extra strict
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 2) then
  set WEIGHTS  = ".002 .005 0.01 0.03 0.05 0.10"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30"
  set JELLY    = "ridg dist sigm 0.10" #Use jelly-body refinement (with looser restraints)
  @ NCYCLE = ($NCYCLE + 5)
  #Use homology restraints
  if (-e $WORKDIR/$PDBID.fasta) then
    set DOHOMOLOGY = 1
  endif
else if ($RESOTYPE == 3) then
  set WEIGHTS  = "0.01 0.03 0.05 0.10 0.30 0.50 0.70"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
else if ($RESOTYPE == 4) then
  set WEIGHTS  = "0.10 0.30 0.50 0.70 1.00 1.50 2.00"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
else if ($RESOTYPE > 4) then
  set WEIGHTS  = "0.70 1.00 1.50 2.00 3.00 5.00"
  set BWEIGHTS = "1.50 1.20 1.00 0.80 0.50 0.30 0.10"
endif

#Set NCS restraints
if ($DONCS == 1) then
  #Check whether NCS is really detected
  if (`grep -c 'Ncsr group' $WORKDIR/${PDBID}_0cyc.log` == 0) then
    #Switch off NCS
    set DONCS      = 0
    set NCSTYPE    =        #Empty NCS command
    set NCSALIGN   =        #Only usefull when NCS is used (alignment cut-off)
    set NCSNEIGH   =        #Only usefull when NCS is used (include neighbouring atoms in restraints)
  endif
endif

#Use harmonic retraints for very small data sets
if ($DOHARMONIC == 1 && $NREFCNT < 1000) then
  set HARMCMD = "ridge atoms 0.05"
  @ NCYCLE = ($NCYCLE + 10)
endif

#Switch off jelly-body refinement
if ($DOJELLY == 0) then
  set JELLY = ""
endif

#Switch off homology based restraints
if ($NOHOMOLOGY == 1) then
  set DOHOMOLOGY = 0
endif

#Use the Wilson B-factor unless it is negative or the resolution is lower than 4.00A; maximise at 50A^2.
setenv BSET  `perl -e 'if ($ENV{URESO} > 3.99 || $ENV{BWILS} < 0) {$bset = 0.5*$ENV{PBAVER}} else {$bset = 0.5*$ENV{BWILS}}; if ($bset > 50) {$bset = 50.00}; if ($bset < 10) {$bset = 10.00}; printf ("%.2f\n", $bset);'`
set TBCMD = "bfac set $BSET"

############################################# Generate external restraints ###############################################
if ( ($GOT_NUC == 'T' && $DOLIBG == 1) || $HBONDREST == 1 || $DOHOMOLOGY == 1 || $DOKRAB == 1) then
  echo "" | tee -a $LOG
  echo "" | tee -a $LOG
  echo "****** Structure specific restraints ******" | tee -a $LOG

  #Run libg three times if the structure has DNA or RNA. PROGRAM: libg
  if ($GOT_NUC == 'T' && $DOLIBG == 1) then
    echo "-Generating nucleic acid restraints" | tee -a $LOG
    libg -p $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/stacking.rest -w sp >  $WORKDIR/libg.log
    libg -p $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/basepair.rest -w bp >> $WORKDIR/libg.log
    libg -p $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/chpucker.rest -w pu >> $WORKDIR/libg.log

    #Make Refmac use command files
    if (-e $WORKDIR/stacking.rest || -e $WORKDIR/basepair.rest || -e $WORKDIR/chpucker.rest) then
      echo "external weight scale 5" > $WORKDIR/nucleic.rest
      if (-e $WORKDIR/stacking.rest) then
        sed 's/  / /g' $WORKDIR/stacking.rest >> $WORKDIR/nucleic.rest
      endif
      if (-e $WORKDIR/basepair.rest) then
        sed 's/  / /g' $WORKDIR/basepair.rest >> $WORKDIR/nucleic.rest
      endif
      if (-e $WORKDIR/chpucker.rest) then
        sed 's/  / /g' $WORKDIR/chpucker.rest >> $WORKDIR/nucleic.rest
      endif  
      echo "external weight scale 1" >> $WORKDIR/nucleic.rest
      set LIBGCMD = "@$WORKDIR/nucleic.rest"
    endif  
  endif

  #Generate a DSSP file
  if ($HBONDREST == 1 || $DOHOMOLOGY == 1) then
    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_0cyc.pdb -o $WORKDIR/${PDBID}_0cyc.dssp >& $WORKDIR/dssp.log

    if ( ! -e $WORKDIR/${PDBID}_0cyc.dssp) then
      #DSSP failed. Cannot make restraints
      echo " o Cannot produce hydrogen bond or homology-based restrains" | tee -a $LOG
      echo "DSSP: general error in restraint generation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                             >> $DEBUG

      #Switch off hydrogen bond and homology restraints
      set HBONDREST  = 0
      set DOHOMOLOGY = 0
  endif

  #Generate hydrogen bond restraints
  if ($HBONDREST == 1) then

    echo "-Generating hydrogen bond restraints" | tee -a $LOG

    #Generate hydrogen bond restraints if a DSSP file can be made.

      #Make the hydrogen bond restraints. PROGRAM: detectHbonds
      $TOOLS/detectHbonds -v \
      -pdb $WORKDIR/${PDBID}_0cyc.pdb \
      -dssp $WORKDIR/${PDBID}_0cyc.dssp \
      -output-name $PDBID \
      -tools $TOOLS > $WORKDIR/hbondrest.log
      if ( ! -e $WORKDIR/${PDBID}_hbonds.rest) then
        echo " o Cannot generate hydrogen bond restraints" | tee -a $LOG
        echo "detectHbonds: general error" >> $DEBUG
        echo "PDB-REDO,$PDBID"             >> $DEBUG
      else
        #Set up restraint commands
        sed 's/ [ ]*/ /g' $WORKDIR/${PDBID}_hbonds.rest > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
        set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
        set HBONDCMD    = "@$WORKDIR/hbond.rest"
      endif
    endif
  endif

  #Generate homology-based restraints
  if ($DOHOMOLOGY == 1 && $GOT_PROT == 'T') then
    echo "-Generating homology-based restraints" | tee -a $LOG

    #Find PDB entries with homologous sequences
    if ($?BLASTP) then
      $BLASTP \
      -evalue 0.001 \
      -num_descriptions 9999 \
      -num_alignments 9999 \
      -query $WORKDIR/$PDBID.fasta \
      -out $WORKDIR/$PDBID.blast \
      -db $TOOLS/pdbredo_seqdb.txt >& $WORKDIR/blast.log
    else
      echo " o BlastP is missing or not configured properly" | tee -a $LOG
      echo " o Cannot generate homology-based restraints"    | tee -a $LOG
      echo "homology: BlastP unavailable" >> $DEBUG
      echo "PDB-REDO,$PDBID"              >> $DEBUG
    endif

    #Make a temporary directory
    mkdir -p $WORKDIR/gbh
    cd $WORKDIR/gbh
    
    #Copy in extra homologous structure models
    if ("$XHOM" != "") then
      echo "-Importing homologous structure models" | tee -a $LOG
      set CNT = 0
      #Loop over files
      foreach FIL ($XHOM)
        #Only use files with exactly 1 CRYST1 card.
        if (`grep -c CRYST1 $FIL` != 1) then
          echo " o Homologous structure model $FIL does not have exactly one CRYST1 card. It will be ignored" | tee -a $LOG
        else
          set CNT  = `expr $CNT + 1`
          set RANK = `seq -w $CNT 99 | head -n 1`
          cp $FIL $WORKDIR/gbh/hom$CNT.pdb
          
          #Construct command line
          set HODERCMD = "$HODERCMD -homol hom$CNT.pdb"
        endif
      end
      echo " o Homologous models imported: $CNT" | tee -a $LOG
    endif

    #Generate the restraints
    $TOOLS/hoder -v \
    -alignout \
    -hitsummary \
    $HODERCMD \
    -pdb $WORKDIR/${PDBID}_0cyc.pdb \
    -dssp $WORKDIR/${PDBID}_0cyc.dssp \
    -blast $WORKDIR/$PDBID.blast \
    -fasta $WORKDIR/$PDBID.fasta \
    -output-name $PDBID \
    -output-dir $WORKDIR/gbh \
    -pdbdir $REDODIR \
    -edsdir $EDSDIR \
    -tools $TOOLS >> homologs.log

    #Set up the restraint commands
    #Homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/gbh/${PDBID}_homologyBased.restr > $WORKDIR/homology.rest #Replace one or more spaces by a single one
    set HOMOLWGTCMD = "EXTERNAL WEIGHT SCALE $HOMOLRESTRWGT"
    set HOMOLCMD    = "@$WORKDIR/homology.rest"

    #Non-homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/gbh/${PDBID}_general.restr > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
    set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
    set HBONDCMD    = "@$WORKDIR/hbond.rest"

    #Go back down
    cd $WORKDIR
  endif
  
  #Make antibody specific restraints
  if ($?KRABSRC && $?PYIGCL && $DOKRAB == 1) then
    echo "-Generating antibody-specific restraints" | tee -a $LOG
    #Source KRAB and create the virtual environment
    source $KRABSRC
      
    #Run KRAB
    krab --pyigclassify_root $PYIGCL $WORKDIR/${PDBID}_0cyc.pdb
    
    #Set up the restraint commands
    if (-e $WORKDIR/krab.rest) then
      #Filter with hbond and homology restraints if they exist
      if (-e $WORKDIR/homology.rest) then
        cp $WORKDIR/homology.rest $WORKDIR/homology.rest.unfiltered
        $TOOLS/filterest.py -v -o $WORKDIR/homology.rest $WORKDIR/krab.rest $WORKDIR/homology.rest.unfiltered > $WORKDIR/filterest.log
      endif
      if (-e $WORKDIR/hbond.rest) then
        cp $WORKDIR/hbond.rest $WORKDIR/hbond.rest.unfiltered
        $TOOLS/filterest.py -v -o $WORKDIR/hbond.rest $WORKDIR/krab.rest $WORKDIR/hbond.rest.unfiltered  >> $WORKDIR/filterest.log
      endif     
      
      set KRABWGTCMD = "EXTERNAL WEIGHT SCALE $KRABRESTRWGT"
      set KRABCMD    = "@$WORKDIR/krab.rest"
    endif
    
    #Append the pepflip skip list
    if (-e $WORKDIR/krab.skip) then
      set BBN_KEEP = "$BBN_KEEP`cat $WORKDIR/krab.skip`"
      set BBN_KEEP = `echo $BBN_KEEP | sed 's/::/:/g'`
    endif
    
    #Switch off the KRAB virtual environment
    deactivate
  endif
endif


############################################### Solvent mask optimisation ################################################


#Swith to simple model
set SOLVENT = SIMP

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Solvent mask optimisation ******" | tee -a $LOG

#Label to go back to after problems with external restraints
solventtest:

#Run Refmac
echo "-Running Refmac grid search" | tee -a $LOG

refmac5 \
XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
XYZOUT $WORKDIR/${PDBID}_solvent.pdb \
HKLIN  $WORKDIR/$PDBID.mtz \
$LIBLIN \
$SCATLIN \
<<eof >$WORKDIR/${PDBID}_solvent.log
  $SCATTERCMD
  make check NONE
  make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
    ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
  refi type REST resi MLKF meth CGMAT bref MIXE
  $REFIRES
  ncyc 0
  scal type $SOLVENT $SCALING
  solvent YES
  solvent optimise
  $NCSSTRICT
  $LOWMEM
  weight $WGTSIG AUTO 2.50
  monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
   chiral 10.0   bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
  $TWIN
  $METALCMD
  $LIBGCMD
  $RESTCMD
  $HOMOLWGTCMD
  $HOMOLCMD
  $HBONDWGTCMD
  $HBONDCMD
  $KRABWGTCMD
  $KRABCMD
  END
eof
if ($status) then
  echo " o Problem with refmac." | tee -a $LOG
  echo "COMMENT: refmac: error solvent optimisation" >> $DEBUG
  echo "PDB-REDO,$PDBID"                             >> $DEBUG

  #Set up fallback values
  set VDWPROBE = 'NA'
  set IONPROBE = 'NA'
  set RSHRINK  = 'NA'

else if (`grep -c 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log` > 0) then
  #There are libg or zen restraint problems report only the first time
  if ($NPRUNE == 0) then
    echo " o Problematic external restraints" | tee -a $LOG
    echo "COMMENT: refmac: external restraint problems" >> $DEBUG
    echo "PDB-REDO,$PDBID"                              >> $DEBUG
  endif

  #try pruning the restraints and try again.
  if ($NPRUNEN < $MAXPRUNE && $NPRUNEM < $MAXPRUNE) then
    cp $WORKDIR/${PDBID}_solvent.log $WORKDIR/${PDBID}_solventv$NPRUNE.log
    if (`grep -c 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log` != 0) then
      set BADREST = `grep -A 1 'Error: At least one of the atoms from the restraints' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | cut -c 6-`
    else if (`grep -c 'Plane number can be 1 or 2' $WORKDIR/${PDBID}_solvent.log` != 0) then
      set BADREST = `grep -A 1 'Plane number can be 1 or 2' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | cut -c 6-`
    else
      echo "   * Cannot locate problematic restraint. Stopping." | tee -a $LOG
      echo "COMMENT: refmac: unsolvable external restraint problems" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
      exit(1)
    endif
    echo "   * Removing problematic restraint" | tee -a $LOG
    if (-e $WORKDIR/nucleic.rest) then
      cp $WORKDIR/nucleic.rest $WORKDIR/temp.rest
      grep -v "$BADREST" temp.rest > $WORKDIR/nucleic.rest
      @ NPRUNEN = ($NPRUNEN + 1)
    endif
    if (-e $WORKDIR/metal.rest) then
      cp $WORKDIR/metal.rest $WORKDIR/temp.rest
      grep -v "$BADREST" temp.rest > $WORKDIR/metal.rest
      @ NPRUNEM = ($NPRUNEM + 1)
    endif
    @ NPRUNE = ($NPRUNEM + $NPRUNEN)
    goto solventtest
  else
    #Give up for the specific restraint type
    if ($NPRUNEM == $MAXPRUNE)  then
      #Stop using metal restraints (reset the counter)
      set METALCMD =
      @ NPRUNEM = ($MAXPRUNE - 1)
    endif
    if ($NPRUNEN == $MAXPRUNE)  then
      #Stop using nucleic acid restraints (reset the counter)
      set LIBGCMD =
      @ NPRUNEN = ($MAXPRUNE - 1)
    endif
    cp $WORKDIR/${PDBID}_solvent.log $WORKDIR/${PDBID}_solventv$NPRUNE.log
    goto solventtest
  endif
else
  #Set up solvent parameters.
  set VDWPROBE = `grep -a -A 3 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $4}' | cut -c 1-3`
  set IONPROBE = `grep -a -A 4 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $4}' | cut -c 1-3`
  set RSHRINK  = `grep -a -A 5 -e 'Minimum R f*ree' $WORKDIR/${PDBID}_solvent.log | tail -n 1 | awk '{print $3}' | cut -c 1-3`
  set MASKPAR  = "solvent vdwprobe $VDWPROBE ionprobe $IONPROBE rshrink $RSHRINK"

  #Report
  echo " o VDW probe: $VDWPROBE" | tee -a $LOG
  echo " o Ion probe: $IONPROBE" | tee -a $LOG
  echo " o Shrinkage: $RSHRINK"  | tee -a $LOG

  #Clean
  rm $WORKDIR/${PDBID}_solvent.pdb
endif


########################################## Extend the resolution #########################################################

#Extend only if the data has higher resolution than the PDB header
if (`perl -e 'if ($ENV{RESO} - $ENV{DRESH} > 0.10) {print "1"} else {print "0"};'` == 1 && $RESOCHECK == 1) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing resolution cut-off ******" | tee -a $LOG

  #Don't come back here
  set RESOCHECK = 0

  #Run binliner
  echo "-Running binliner" | tee -a $LOG
  #PROGRAM: binliner
  $TOOLS/binliner -v \
  $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA \
  $WORKDIR/$PDBID.cif \
  $RESO $DATARESH > $WORKDIR/binliner.log
  echo "-Testing resolution cut-offs:" `tail -n 1 $WORKDIR/binliner.log` | tee -a $LOG

  #Initialise values
  set RLOWER   = 1.00
  set RFLOWER  = 1.00
  set WRFLOWER = 1.00
  set FLLOWER  = 999999.9
  set FCCLOWER = 0.00
  set RESLOWER = $RESO
  set RESSTEP  = -1
  @ NRCYCLE = ($NCYCLE + 10)

  #Start a pretty logfile
  echo "CUTOFF RLOWER RHIGHER RFLOWER RFHIGHER WRFLOWER WRFHIGHER FLLOWER FLHIGHER FCCLOWER FCCHIGHER" > $WORKDIR/${PDBID}_resotest.log

  #Refine to extend the resolution
  foreach RESCO (`tail -n 1 $WORKDIR/binliner.log`)

    echo " o Testing resolution $RESCO" | tee -a $LOG
    #Set resolution cut-offs
    set REFIRES = "RESO $RESCO"

    #Refine against data
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
    XYZOUT $WORKDIR/${PDBID}_res$RESCO.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_all$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      $TBCMD
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc $NRCYCLE
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $LIBGCMD
      $RESTCMD
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      $KRABWGTCMD
      $KRABCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Calculate R(-free) against lower resolution
    set TESTRES = "RESO $RESLOWER"

    #Run Refmac 0-cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_res$RESCO.pdb \
    XYZOUT $WORKDIR/${PDBID}_restest_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest_0cyc.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_high2low$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen YES hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref ISOT
      $TESTRES
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $LIBGCMD
      $RESTCMD
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      $KRABWGTCMD
      $KRABCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Accept or reject the higher resolution data
    set RHIGHER   = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $2}'`
    set RFHIGHER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $3}'`
    set WRFHIGHER = `grep "Free weighted R2 factor" $WORKDIR/${PDBID}_high2low$RESCO.log | awk '{print $6}'`
    set FLHIGHER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_high2low$RESCO.log | head -n 1 | awk '{print $6}'`
    set FCCHIGHER = `grep "Free correlation coefficient" $WORKDIR/${PDBID}_high2low$RESCO.log | awk '{print $5}'`

    #Continue logfile
    echo $RESCO $RLOWER $RHIGHER $RFLOWER $RFHIGHER $WRFLOWER $WRFHIGHER $FLLOWER $FLHIGHER $FCCLOWER $FCCHIGHER >> $WORKDIR/${PDBID}_resotest.log

    #Is the best resolution so far the same as the higher resolution? PROGRAM: resolute
    if (`$TOOLS/resolute $WORKDIR/${PDBID}_resotest.log 1` == $RESCO) then
      #Yes. Current resolution is better than previous one. Continue.
      @ RESSTEP = ($RESSTEP + 1)
    else
      #No. Reject the resolution and exit the loop.
      goto gotbestreso
    endif

    #Prepare for next cycle get R(-free) from 0 cycles
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_res$RESCO.pdb \
    XYZOUT $WORKDIR/${PDBID}_restest_0cyc.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_restest_0cyc.mtz \
    $LIBLIN \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_0cyc$RESCO.log
      $SCATTERCMD
      make check NONE
      make hydrogen YES hout YES peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref ISOT
      $REFIRES
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG AUTO 2.50
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $HARMCMD
      $METALCMD
      $LIBGCMD
      $RESTCMD
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      $KRABWGTCMD
      $KRABCMD
      END
eof
    if ($status) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in resolution extension" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif

    #Set values for next round
    set RLOWER   = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $2}'`
    set RFLOWER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $3}'`
    set WRFLOWER = `grep "Free weighted R2 factor" $WORKDIR/${PDBID}_0cyc$RESCO.log | awk '{print $6}'`
    set FLLOWER  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_0cyc$RESCO.log | head -n 1 | awk '{print $6}'`
    set FCCLOWER = `grep "Free correlation coefficient" $WORKDIR/${PDBID}_0cyc$RESCO.log | awk '{print $5}'`
    setenv URESO $RESCO
    set RESLOWER = $RESCO
  end

  #Update the reflection count
gotbestreso:

  #Run MTZUTILS to cut the data
  mtzutils \
  HKLIN  $WORKDIR/$PDBID.mtz \
  HKLOUT $WORKDIR/lowres.mtz \
  <<eof >>$WORKDIR/mtz_creation.log
    RESOLUTION $URESO $DATARESL
    END
eof
  if ($status) then
    echo "   * Error using MTZUTILS. Cannot continue." | tee -a $LOG
    echo "COMMENT: mtzutils: general error (2)" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                      >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    exit(1)
  endif

  #Get the correct number of (R-free) reflections
  mtz2various \
  HKLIN $WORKDIR/lowres.mtz \
  HKLOUT $WORKDIR/temp.hkl \
<<eof >>$WORKDIR/mtz_creation.log
    OUTPUT CIF data_temp
    LABIN FP=FP SIGFP=SIGFP FREE=FREE
    END
eof
  set NTSTCNT = `grep -cE ' f ' $WORKDIR/temp.hkl`
  setenv PNTSTCNT $NTSTCNT  #Set value for bits of Perl code in the script
  set NREFCNT = `grep -cE ' [of] ' $WORKDIR/temp.hkl`
  setenv PWORKCNT `grep -cE ' o ' $WORKDIR/temp.hkl`

  #Cleanup
  rm $WORKDIR/temp.hkl
  rm $WORKDIR/lowres.mtz
  rm $WORKDIR/${PDBID}_restest.mtz
  rm $WORKDIR/${PDBID}_restest_0cyc.pdb
  rm $WORKDIR/${PDBID}_restest_0cyc.mtz

  #Update the resolution based settings
  echo "-High resolution cut-off: $URESO" | tee -a $LOG
  set REFIRES = "RESO $URESO"

  #Update the number of refinement cycles:
  if ($RESSTEP > 0) then
    set TCYCLE = $NCYCLE
    @ NCYCLE = ($NCYCLE + $RESSTEP * 5)
    echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  endif

  #Update R-factor stats
  set RCAL  = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_all$URESO.log | tail -n 1 | awk '{print $2}'`
  set RFCAL = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_all$URESO.log | tail -n 1 | awk '{print $3}'`
  setenv PRCAL  $RCAL
  setenv PRFCAL $RFCAL

  goto resobasedsettings
endif


####################################### Decide on type of B-factor refinement ############################################

#Try refinement with anisotropic B's and no TLS if the data parameter ratio allows it (not with heavy strict NCS)
if ($STRICTNCS == 1 && $PMTRIX > 9) then
  set BREFTYPE     = "ISOT"
  set BMODELPHRASE = 'only an isotropic B-factor model was considered.'
else
  set BREFTYPE = `perl -e 'if ($ENV{URESO} > 1.94) {print "ISOT\n"} elsif ($ENV{PWORKCNT}/$ENV{PATMCNT} > 30.0) {print "ANISOT\n"} elsif ($ENV{PWORKCNT}/$ENV{PATMCNT} > 13.0) {print "TEST\n"} else {print "ISOT\n"};'`
endif
set REFPATM  = `perl -e 'printf ("%.1f\n", $ENV{PWORKCNT}/$ENV{PATMCNT});'`

#Test the best B-factor restraint type
if ($BREFTYPE == TEST) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing B-factor refinement type ******" | tee -a $LOG

  #Run refmac with and without anisotropic B-factors.
  foreach TYPETEST (`echo "ANISOT ISOT"`)

    #Return label for job launching
anisooriso:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log
    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      echo "-Testing ${TYPETEST}ropic B-factors" | tee -a $LOG

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
      XYZOUT $WORKDIR/${PDBID}_${TYPETEST}ropic.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz \
      $LIBLIN \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${TYPETEST}ropic.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        $TBCMD
        refi type REST resi MLKF meth CGMAT bref $TYPETEST
        $REFIRES
        ncyc 50
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bspherem 10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $LIBGCMD
        $RESTCMD
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy remarks 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KRABWGTCMD
        $KRABCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto anisooriso
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists; then clean up
  foreach TYPETEST (`echo "ANISOT ISOT"`)
    if (! -e $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz) then
      echo "o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in first B-factor type selection" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                         >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      rm $WORKDIR/${PDBID}_${TYPETEST}ropic.mtz
      rm $WORKDIR/${PDBID}_${TYPETEST}ropic.pdb
    endif
  end

  #Compare refinement results and decide on the B-factor type. PROGRAM: bselect
  set BREFTYPE = `$TOOLS/bselect $WORKDIR/${PDBID}_ANISOTropic.log $WORKDIR/${PDBID}_ISOTropic.log`

endif

#Setup the refinement
if ($BREFTYPE == "ISOT") then
  #Set number of parameters per atom
  setenv PPPATM  `echo "4"`

  #Use TLS unless it is surpressed
  if ($DOTLS == 1) then

    #Set command file
    if(`ls ????.tls | wc -l` != 0) then
      set TLSCMD  = `echo refi tlsc $TLSCYCLE`

      #Keep only the TLS group definitions, recalculate the origin and the tensors
      if (-e $WORKDIR/${PDBID}.tls) then
        cp $WORKDIR/${PDBID}.tls $WORKDIR/${PDBID}.tls_original
        grep -E 'TLS|RANGE|^ ' $WORKDIR/${PDBID}.tls_original > $WORKDIR/${PDBID}.tls
      endif
    endif
  endif

  #Create B-factor model phrase for the summary
  set BMODELPHRASE = 'both an isotropic and an anisotropic B-factor model were considered, and the isotropic B-factor model was selected based on the <a href="http://www.cmbi.ru.nl/pdb_redo/faq.html#hamr" target="_blank">Hamilton R ratio test</a>.'
else
  #Set number of parameters per atom
  setenv PPPATM  `echo "9"`

  #Do not use TLS
  set DOTLS    = 0  #No TLS at all
  cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb

  #Do extra cycles in B-restraint weight optimisation
  set BTESTCYC = 15

  #Create B-factor model phrase for the summary
  set BMODELPHRASE = 'both an isotropic and an anisotropic B-factor model were considered, and the anisotropic B-factor model was selected based on the <a href="http://www.cmbi.ru.nl/pdb_redo/faq.html#hamr" target="_blank">Hamilton R ratio test</a>.'
endif


################################################### Show refinement setting ##############################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Refinement settings ******" | tee -a $LOG

#Type of refinement
echo "-B-factor model" | tee -a $LOG
echo " o Number of atoms      : $PATMCNT"  | tee -a $LOG
echo " o Number of reflections: $PWORKCNT" | tee -a $LOG
echo " o Reflections per atom : $REFPATM"  | tee -a $LOG
echo " o B-factor type        : ${BREFTYPE}ropic" | tee -a $LOG
if ($DOTLS == 0) then
  echo " o TLS-models           : not used" | tee -a $LOG
else
  echo " o TLS-models           : used" | tee -a $LOG
endif

#NCS restraint type
echo "-Non-crystallographic symmetry" | tee -a $LOG
if ($DONCS == 1 && $STRICTNCS == 1) then
  echo " o Using local NCS restraints and global NCS constraints" | tee -a $LOG
else if ($DONCS == 1) then
  echo " o Using local NCS restraints" | tee -a $LOG
else if ($STRICTNCS == 1) then
  echo " o Using NCS constraints" | tee -a $LOG
else
  echo " o No NCS restraints/constraints used" | tee -a $LOG
endif

#Twinning
echo "-Twinning" | tee -a $LOG
if ($DOTWIN == 1) then
  if ($TWIN == 'twin') then
    echo " o Detwinning during refinement" | tee -a $LOG
  else
    echo " o No twinning detected" | tee -a $LOG
  endif
else
  echo " o Detwinning surpressed"  | tee -a $LOG
endif

#No R-free was reported so a new target must be defined. The 'Unbiased' R-free is used.
if ($RFHEAD == 0) then
  echo "-No R-free reported in PDB header" | tee -a $LOG
  echo " o R-free target is now $RFCALUNB" | tee -a $LOG
endif

#Warn for high structure RMSZ scores
if ($RMSZB == 'huge' || $RMSZA == 'huge' || `perl -e 'if ($ENV{PRMSZB} > 1.000) {print "1"} elsif ($ENV{PRMSZA} > 1.000) {print "1"} else {print "0"};'` == 1) then
  echo "-High Refmac bond length or bond angle RMSZ detected" | tee -a $LOG
  echo " o Bond length RMSZ: $RMSZB"                          | tee -a $LOG
  echo " o Bond angle RMSZ : $RMSZA"                          | tee -a $LOG
  echo " o Geometric targets will be relaxed"                 | tee -a $LOG
endif

#Check for strangely low R-free values (only when the original R-free set was used)
if ($GOTR == 1) then
  #Check for high Z-score or R-free lower than a certain cut-off (0.0 or R+0.33(Rfree-R))
  if ($RHEAD == 1 && $RFHEAD == 1) then
    set ZCALERR = `perl -e 'if ($ENV{PRFCALZ} > 7.0) {print "1\n"} elsif (($ENV{PRFCAL} - $ENV{PRCAL}) < 0.33*($ENV{PRFREE} - $ENV{PRFACT})) {print "1\n"} else {print "0\n"};'`
  else
    set ZCALERR = `perl -e 'if ($ENV{PRFCALZ} > 7.0) {print "1\n"} elsif (($ENV{PRFCAL} - $ENV{PRCAL}) < 0.0) {print "1\n"} else {print "0\n"};'`
  endif
  if ($ZCALERR == 1) then
    echo "-Severe R-free bias and possible test set problem!" | tee -a $LOG
    #Do 10 cycles extra refinement
    set TCYCLE = $NCYCLE
    @ NCYCLE = ($NCYCLE + 10)
    echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  endif
else
  set ZCALERR = 0
endif

#Compensate for new R-free set
if ($GOTR == 0) then
  echo "-New R-free set" | tee -a $LOG
  echo " o Calculated R and 'Unbiased' R-free values will be used for reference" | tee -a $LOG

  #Use 0.5*Wilson B-factor unless the resolution lower than 4.00A; maximise at 50A^2.
  set BCMD = "bfac set $BSET"
  echo " o Resetting B-factors to $BSET to remove model bias" | tee -a $LOG

  #Do 10 restrained refinement cycles extra
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 10)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
endif

#Compensate for legacy model
if ($LEGACY == 1 ) then
  echo "-Legacy mode" | tee -a $LOG
  echo " o Using calculated R-factor ($RCAL) as refinement target" | tee -a $LOG

  #Do 10 restrained refinement cycles extra
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 10)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
endif

#Compensate for anisotropic B-factors
if ($BREFTYPE == "ANISOT") then
  echo "-Anisotropic B-factor refinement" | tee -a $LOG
  #Do extra cycles of refinement
  set TCYCLE = $NCYCLE
  @ NCYCLE = ($NCYCLE + 20)
  echo " o Doing $NCYCLE cycles of restrained refinement instead of $TCYCLE" | tee -a $LOG
  echo " o Updating R-factor and R-free details"
endif

#Warn about occupancy refinement
if ($OCCREF > 0) then
  echo "-Occupancy refinement" | tee -a $LOG
  echo " o Performing occupancy refinement on $OCCREF residues" | tee -a $LOG
endif

#warn about external restraints
if ($METALCMD != "" || $LIBGCMD != "" || $RESTCMD != "" || $HBONDCMD != "" || $KRABCMD != "" || $HOMOLCMD != "" || `echo $JELLY | grep -c ridg` != 0 || `echo $HARMCMD | grep -c ridg` != 0) then
  echo "-Additional restraints" | tee -a $LOG
  if (`echo $JELLY | grep -c ridg` != 0) then
    echo " o Using jelly-body restraints" | tee -a $LOG
  endif
  if (`echo $HARMCMD | grep -c ridg` != 0) then
    echo " o Using harmonic restraints" | tee -a $LOG
  endif
  if (-e $WORKDIR/metal.rest) then
    set NMETALREST = `grep -c exte $WORKDIR/metal.rest`
    if ($NMETALREST > 0) then
      echo " o Using $NMETALREST metal site restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/nucleic.rest) then
    set NNUCLEICREST = `grep -c exte $WORKDIR/nucleic.rest`
    if ($NNUCLEICREST > 0) then
      echo " o Using $NNUCLEICREST DNA/RNA restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/homology.rest) then
    set NHOMOLREST = `grep -c exte $WORKDIR/homology.rest`
    if ($NHOMOLREST > 0) then
      echo " o Using $NHOMOLREST homology restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/hbond.rest) then
    set NHBONDREST = `grep -c exte $WORKDIR/hbond.rest`
    if ($NHBONDREST > 0) then
      echo " o Using $NHBONDREST hydrogen bond restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/krab.rest) then
    set NKRABREST = `grep -c exte $WORKDIR/krab.rest`
    if ($NKRABREST > 0) then
      echo " o Using $NKRABREST antibody-specific restraints" | tee -a $LOG
    endif
  endif
  if (-e $WORKDIR/external.rest) then
    set NEXTERNALREST = `grep -c exte $WORKDIR/external.rest`
    #Subtract 2 instoduced statements
    @ NEXTERNALREST = ($NEXTERNALREST - 2)
    if ($NEXTERNALREST > 0) then
      echo " o Using $NEXTERNALREST user-provided restraints" | tee -a $LOG
    endif
  endif
endif

#Update the R(-free) statistics if the data/parameter ratio changed
if ($BREFTYPE == "ANISOT" || $URESO != $RESO) then
  #Calculate sigma(R-free)
  set SIGRFCAL = `perl -e 'printf ("%.4f\n", $ENV{PRFCAL}/sqrt($ENV{PNTSTCNT}));'`
  setenv PSIGRFCAL $SIGRFCAL

  #Calculate unbiased R-free/R ratio (thanks to Ian Tickle)
  set RFRRAT = `perl -e '$x = $ENV{PATMCNT}/$ENV{PWORKCNT}; $z = 1-$ENV{PPPATM}*$x/(1+2.5*$x); $z = 0 if $z<0; $a = $ENV{PPPATM}-2.5*(1-$z**5); $x = $a*$x; $rfree_rat = 1-$x; if($rfree_rat <= 0) {$rfree_rat = 1.010} else {$rfree_rat = sqrt((1+$x)/$rfree_rat)} $rfree_rat = 1.200000 if ($ENV{URESO}>2.65 and $rfree_rat>1.200000); $rfree_rat = 1.200000 if ($ENV{URESO}>3.0 and $rfree_rat<1.011); $rfree_rat = 1.450000 if $rfree_rat>1.454; printf ("%.6f\n", $rfree_rat);'`
  setenv PRFRRAT $RFRRAT

  #Calculate unbiased R-free
  set RFCALUNB = `perl -e 'printf ("%.4f\n", $ENV{PRCAL}*$ENV{PRFRRAT});'`
  setenv PRFCALUNB $RFCALUNB

  #R-free Z-score (comparison of R-free with its expected value)
  set RFCALZ = `perl -e 'printf ("%.2f\n", ($ENV{PRFCALUNB}-$ENV{PRFCAL})/$ENV{PSIGRFCAL});'`
  setenv PRFCALZ $RFCALZ

  #Print values
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** R-factor and R-free details (updated) ******" | tee -a $LOG
  echo "Calculated R     : $RCAL"     | tee -a $LOG
  echo "Calculated R-free: $RFCAL"    | tee -a $LOG
  echo "'Unbiased' R-free: $RFCALUNB" | tee -a $LOG
  echo "sigma(R-free)    : $SIGRFCAL" | tee -a $LOG
  echo "R-free Z-score   : $RFCALZ"   | tee -a $LOG
endif

################################################ TLS group optimisation ##################################################

#Only do this if TLS refinement isn't surpressed
if ($DOTLS == 1) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** TLS group optimisation ******" | tee -a $LOG

  #Optimise TLS groups iff any definitions exist.
  set NTLS = `find $WORKDIR -name "*.tls" | wc -l`
  if ($NTLS != 0) then
    echo "-Testing $NTLS TLS group definition(s)" | tee -a $LOG

    #Maximise the used CPU time
    limit cputime 24h

    #Run refmac with different TLS group configurations
    foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)

      #Return label for job launching
ttestrunning:

      #Only launch new jobs when the number of cores is not exceeded
      #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
      jobs > $WORKDIR/jobs.log
      if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

        #Run refmac
        refmac5 \
        XYZIN  $WORKDIR/${PDBID}_0cyc.pdb \
        XYZOUT $WORKDIR/${PDBID}_ttest$TLSG.pdb \
        HKLIN  $WORKDIR/$PDBID.mtz \
        HKLOUT $WORKDIR/${PDBID}_ttest$TLSG.mtz \
        $LIBLIN \
        TLSIN $WORKDIR/$TLSG.tls TLSOUT $WORKDIR/${TLSG}_out.tls \
        $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_ttest$TLSG.log &
          $SCATTERCMD
          make check NONE
          make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
           ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
          $TBCMD
          refi type REST resi MLKF meth CGMAT bref ISOT
          $REFIRES
          $TLSCMD
          tlsd waters exclude
          ncyc 0
          scal type $SOLVENT $SCALING
          solvent YES
          $MASKPAR
          $LOWMEM
          weight AUTO
          monitor MEDIUM -
            torsion 10.0 distance 10.0 angle 10.0 plane 10.0 chiral 10.0 -
            bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
          $NCSTYPE
          $NCSALIGN
          $NCSNEIGH
          $NCSSTRICT
          $TWIN
          blim 2.0 999.0
          labin  FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
          $ANOMCMD
          pdbout copy remarks 280 350
          NOHARVEST
          $HOMOLWGTCMD
          $HOMOLCMD
          $HBONDWGTCMD
          $HBONDCMD
          $KRABWGTCMD
          $KRABCMD
          kill $TOOLS/pdb_redo.refmac
          END
eof
      else
        #Wait a bit to start again
        sleep 10
        goto ttestrunning
      endif
    end

    #Wait for the jobs to finish
    wait

    #Unset the CPU time limit
    limit cputime unlimited

    #Check for errors and report results
    foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
      if (! -e $WORKDIR/${TLSG}_out.tls) then
        echo " o Problem with Refmac using $TLSG.tls" | tee -a $LOG
        mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
        @ NTLS = $NTLS - 1 
      else
        #Mine out the R-free and report
        set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1 | awk '{print $3}'`
        echo " o Tested tls groups from $TLSG.tls (R-free = $TFREE)" | tee -a $LOG
      endif
    end

    #Remove TLS group definitions that cause crazy values in the tensors or B-factors
    if ($NTLS != 0) then
      echo "-Filtering results" | tee -a $LOG
      foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
        if (`grep -c -E '\*{6}' ${TLSG}_out.tls` != 0) then
          echo " o Problem with TLS group definition $TLSG.tls" | tee -a $LOG
          mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
          @ NTLS = $NTLS - 1
        endif
        if (-e $WORKDIR/$TLSG.tls && `grep '^[AH][TE]' ${PDBID}_ttest${TLSG}.pdb | grep -c -E '\*{6}'` != 0) then
          echo " o Problem with TLS group definition $TLSG.tls" | tee -a $LOG
          mv $WORKDIR/$TLSG.tls $WORKDIR/$TLSG.notls
          @ NTLS = $NTLS - 1
        endif
      end
    endif  

    #Create .ttest file for picker iff any TLS groups are approved
    if ($NTLS != 0) then
      #Grep out the refinement stats, and take the first data line to obtain R(-free) targets. Use the first valid log file.
      set VALID1 = `ls -Sr ????.tls | head -n 1 | cut -c 1-4`
      set TTR  = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_ttest$VALID1.log | tail -n 1 | awk '{print $2}'`
      set TTRF = `grep -a -A 2 'Ncyc    Rfact    Rfree     FOM' $WORKDIR/${PDBID}_ttest$VALID1.log | tail -n 1 | awk '{print $3}'`

      echo "$TTR $TTRF" > $WORKDIR/${PDBID}.ttest   #Use R and R-free obtained after resetting the B-factor
      foreach TLSG (`ls -Sr ????.tls | cut -c 1-4`)
        if (-e $WORKDIR/${PDBID}_ttest$TLSG.pdb) then
          if ($TLSG == $VALID1) then
            set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1`
            echo "$TLSG $LINE" >> $WORKDIR/${PDBID}.ttest
          else
            #Do a hamilton test first
            if (`$TOOLS/bselect -t $WORKDIR/${PDBID}_ttest$VALID1.log $WORKDIR/${PDBID}_ttest$TLSG.log` != 'LOG1') then
              set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_ttest$TLSG.log | head -n 1`
              echo "$TLSG $LINE" >> $WORKDIR/${PDBID}.ttest
            else
              echo " o TLS group definition $TLSG.tls causes over-fitting" | tee -a $LOG
            endif
          endif
        else
          echo " o $WORKDIR/${PDBID}_ttest$TLSG.pdb was missing" | tee -a $LOG
          set TTEST = 1
        endif
      end

      #Pick the best TLS group configuration. PROGRAM: picker
      set OPTTLSG = `$TOOLS/picker $WORKDIR/${PDBID}.ttest $NTSTCNT $RFRRAT 10.0 10.0` #The geometry is ignored here

      #Also pick second best TLS group.
      if (`ls -Sr ????.tls | wc -l` != 1) then
        grep -v $OPTTLSG $WORKDIR/${PDBID}.ttest > $WORKDIR/${PDBID}.ttest2
        set OPTTLSG2 = `$TOOLS/picker -s $WORKDIR/${PDBID}.ttest2 $NTSTCNT $RFRRAT 10.0 10.0`
      else
        set OPTTLSG2 = 'none'
      endif

      #Print values
      if ($OPTTLSG == 'none') then
        echo "-TLS does not seem to work for this structure" | tee -a $LOG
        echo " o TLS refinement will not be used"  | tee -a $LOG
        #Set up input for B-weight optimisation
        cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
        set DOTLS    = 0
        set TLSCMD   =        #Empty TLS command
        set TLSFILS  =        #No TLS files specified
      else
        echo "-TLS groups in $OPTTLSG.tls will be used in refinement" | tee -a $LOG
        #Set up input for B-weight optimisation
        cp $WORKDIR/${OPTTLSG}_out.tls         $WORKDIR/${PDBID}_refin.tls
        set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`
        set BCMD    =  #No need to reset the B-factors anymore
        #Run TLSanl to get a proper input PDB file if an older Refmac was used
        tlsanl \
        XYZIN $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb \
        XYZOUT $WORKDIR/${PDBID}_refin.pdb \
<<eof >& $WORKDIR/${PDBID}_tlsanl.log
        BINPUT t
        BRESID f
        ISOOUT RESI
        NUMERICAL
        END
eof
        if($status) then
          echo " o Problem with TLSanl for TLS groups in $OPTTLSG.tls" | tee -a $LOG
          echo "COMMENT: TLSanl: general error for optimal TLS group" >> $DEBUG
          echo "PDB-REDO,$PDBID"                                      >> $DEBUG
          #Try using the second best TLS group configuration
          if ($OPTTLSG2 != 'none') then

            set OPTTLSG = $OPTTLSG2
            echo " o TLS groups in $OPTTLSG.tls will be used in refinement instead" | tee -a $LOG

            #Run TLSanl again to get a proper input PDB file
            tlsanl \
            XYZIN $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb \
            XYZOUT $WORKDIR/${PDBID}_refin.pdb \
<<eof >$WORKDIR/${PDBID}_tlsanl2.log
            BINPUT t
            BRESID f
            ISOOUT RESI
            NUMERICAL
            END
eof
            if($status) then #second failure
              echo " o Problem with TLSanl for TLS groups in $OPTTLSG.tls" | tee -a $LOG
              echo "COMMENT: TLSanl: general error for second best TLS group" >> $DEBUG
              echo "PDB-REDO,$PDBID"                                          >> $DEBUG
              cp $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb $WORKDIR/${PDBID}_refin.pdb
              #The output files will have the total B-factor instead of the residual. Circumvent by resetting the B-factors.
              set BCMD = `echo $TBCMD`
            endif

            #Set up input for B-weight optimisation (again)
            cp $WORKDIR/${OPTTLSG}_out.tls         $WORKDIR/${PDBID}_refin.tls
            set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`

          else
            rm $WORKDIR/${PDBID}.ttest2
            cp $WORKDIR/${PDBID}_ttest$OPTTLSG.pdb $WORKDIR/${PDBID}_refin.pdb
            #The output files will have the total B-factor instead of the residual. Circumvent by resetting the B-factors.
            set BCMD = `echo $TBCMD`
          endif
        endif
      endif
    else
      #All TLS definitions failed
      echo "-No usable TLS group definitions"   | tee -a $LOG
      echo " o TLS refinement will not be used" | tee -a $LOG
      set DOTLS   = 0  #Do not use TLS
      set TLSCMD  =    #Empty TLS command
      set TLSFILS =    #No TLS files specified
      cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
    endif
  else
    #No TLS groups could be made
    echo "-Cannot create any TLS group definitions" | tee -a $LOG
    echo " o TLS refinement will not be used"       | tee -a $LOG
    echo "COMMENT: Could not create any TLS groups" >> $DEBUG
    echo "PDB-REDO,$PDBID"                          >> $DEBUG
    set DOTLS   = 0  #Do not use TLS
    set TLSCMD  =    #Empty TLS command
    set TLSFILS =    #No TLS files specified
    cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
  endif

else
  #Not using TLS ensure that the input file for the B-weight optimisation exists
  cp $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/${PDBID}_refin.pdb
endif

##################################### Individual B-factors or one overall B-factor? ######################################

#Try refinement with one overall B if TLS works and there are less than 4 reflections per atom.
if ($STRICTNCS == 1 && $PMTRIX > 9) then
  set BTYPE = "ISOT"
  set BMODELPHRASE = 'only an isotropic B-factor model was considered.'
else
  set BTYPE = `perl -e 'if ($ENV{PWORKCNT}/$ENV{PATMCNT} < 4) {print "TEST\n"} else {print "ISOT\n"};'`
endif

#Test the best B-factor restraint type
if ($BTYPE == TEST) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Testing B-factor refinement type ******" | tee -a $LOG

  #Set the B-factor restraint weight
  set TIGHTB = `echo $BWEIGHTS | cut -c 1-4`

  #Set all B-factors to a single value if no TLS is used
  if ($DOTLS == 0) then
    set BBCMD = `echo $TBCMD`
  else
    set BBCMD =  #No B-factor resetting (the B-factors were reset during the TLS optimisation)
  endif

  #Run refmac with and without individual B-factors.
  foreach TYPETEST (`echo "ISOT OVER"`)

    #Return label for job launching
isoorover:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log
    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      #Report the type of refinement
      if ($TYPETEST == "ISOT") then
        echo "-Testing isotropic B-factors"  | tee -a $LOG
      else
        echo "-Testing one overall B-factor" | tee -a $LOG
      endif

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_refin.pdb \
      XYZOUT $WORKDIR/${PDBID}_${TYPETEST}.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_${TYPETEST}.mtz \
      $LIBLIN \
      $TLSFILS \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${TYPETEST}.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        $BBCMD
        refi type REST resi MLKF meth CGMAT bref $TYPETEST
        $REFIRES
        ncyc $NCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $LIBGCMD
        $RESTCMD
        temp $TIGHTB
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy remarks 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KRABWGTCMD
        $KRABCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto isoorover
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists; then clean up
  foreach TYPETEST (`echo "ISOT OVER"`)
    if (! -e $WORKDIR/${PDBID}_${TYPETEST}.mtz) then
      echo "o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in second B-factor type selection" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                          >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      rm $WORKDIR/${PDBID}_${TYPETEST}.mtz
      rm $WORKDIR/${PDBID}_${TYPETEST}.pdb
    endif
  end

  #Compare refinement results and decide on the B-factor type
  set BREFTYPE = `$TOOLS/bselect $WORKDIR/${PDBID}_ISOT.log $WORKDIR/${PDBID}_OVER.log`

  #Setup the refinement
  if ($BREFTYPE == "OVER") then
    echo "-PDB-REDO will use one overall B-factor in refinement" | tee -a $LOG

    #Set number of parameters per atom
    setenv PPPATM  `echo "3"`

    #Set all B-factors to a single value if no TLS is used
    if ($DOTLS == 0) then
      set BCMD = `echo $TBCMD`
    endif

    #Re-calculate unbiased R-free/R ratio (thanks to Ian Tickle)
    set RFRRAT = `perl -e '$x = $ENV{PATMCNT}/$ENV{PWORKCNT}; $z = 1-$ENV{PPPATM}*$x/(1+2.5*$x); $z = 0 if $z<0; $a = $ENV{PPPATM}-2.5*(1-$z**5); $x = $a*$x; $rfree_rat = 1-$x; if($rfree_rat <= 0) {$rfree_rat = 1.010} else {$rfree_rat = sqrt((1+$x)/$rfree_rat)} $rfree_rat = 1.200000 if ($ENV{URESO}>2.65 and $rfree_rat>1.200000); $rfree_rat = 1.200000 if ($ENV{URESO}>3.0 and $rfree_rat<1.011); $rfree_rat = 1.450000 if $rfree_rat > 1.454; printf ("%.6f\n", $rfree_rat);'`
    setenv PRFRRAT $RFRRAT

    #Calculate unbiased R-free
    set RFCALUNB = `perl -e 'printf ("%.4f\n", $ENV{PRCAL}*$ENV{PRFRRAT});'`
    setenv PRFCALUNB $RFCALUNB

    #R-free Z-score (comparison of R-free with its expected value)
    set RFCALZ = `perl -e 'printf ("%.2f\n", ($ENV{PRFCALUNB}-$ENV{PRFCAL})/$ENV{PSIGRFCAL});'`
    setenv PRFCALZ $RFCALZ

    #Create B-factor model phrase for the summary
    set BMODELPHRASE = 'both an isotropic and a flat model with one overall B-factor were considered, and the flat B-factor model was selected based on the <a href="http://www.cmbi.ru.nl/pdb_redo/faq.html#hamr" target="_blank">Hamilton R ratio test</a>.'
  else
    echo "-PDB-REDO will use isotropic B-factors in refinement" | tee -a $LOG

    #Create B-factor model phrase for the summary
    set BMODELPHRASE = 'both an isotropic and a flat model with one overall B-factor were considered, and the isotropic B-factor model was selected based on the <a href="http://www.cmbi.ru.nl/pdb_redo/faq.html#hamr" target="_blank">Hamilton R ratio test</a>.'
  endif
endif

################################################# B-weight optimisation ##################################################

if ($BREFTYPE == ISOT || $BREFTYPE == ANISOT) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** B-weight optimisation ******" | tee -a $LOG
  echo "-Testing B-factor restraint weights: $BWEIGHTS" | tee -a $LOG

  #Run refmac with predefined B-factor restraint weights
  foreach BWGT (`echo $BWEIGHTS`)

    #Return label for job launching
bwgtrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      #Run Refmac. Errors are caught later by checking the output.
      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_refin.pdb \
      XYZOUT $WORKDIR/${PDBID}_btest$BWGT.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_btest$BWGT.mtz \
      $LIBLIN \
      $TLSFILS \
      $SCATLIN \
<<eof >$WORKDIR/${PDBID}_btest$BWGT.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        $BCMD
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        ncyc $BTESTCYC
        tlsd waters exclude
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG AUTO 2.50
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $LIBGCMD
        $RESTCMD
        temp $BWGT
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy remarks 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KRABWGTCMD
        $KRABCMD
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto bwgtrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists
  foreach BWGT (`echo $BWEIGHTS`)
    if (! -e $WORKDIR/${PDBID}_btest$BWGT.mtz) then
      echo " o Problem with refmac. Cannot continue." | tee -a $LOG
      echo "COMMENT: refmac: error in B-weight optimisation" >> $WHYNOT
      echo "PDB-REDO,$PDBID"                                 >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    endif
  end

  #Report
  foreach BWGT (`echo $BWEIGHTS`)
    set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_btest$BWGT.log | head -n 1 | awk '{print $3}'`
    echo " o Performed B-weight test with weight $BWGT (R-free = $TFREE)" | tee -a $LOG
  end

  #Create .btest file for picker
  echo "-Selecting best B weight" | tee -a $LOG
  set BTEST = 0
  if ($GOTR == 0 || $ZCALERR == 1) then
    echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}.btest  #Use recalculated R and 'unbiased' R-free as benchmarks
  else
    echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}.btest  #Use R and R-free obtained from the recalculation as benchmarks
  endif
  foreach BWGT (`echo $BWEIGHTS`)
    if (-e $WORKDIR/${PDBID}_btest$BWGT.pdb) then
      set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_btest$BWGT.log | head -n 1`
      echo "$BWGT $LINE" >> $WORKDIR/${PDBID}.btest
    else
      echo " o $WORKDIR/${PDBID}_btest$BWGT.pdb was missing." | tee -a $LOG
      set BTEST = 1
    endif
  end

  #Pick the best weight. If no weight is found 'none' is returned and the re-refinement runs with default settings.
  set BBEST = `$TOOLS/picker -s $ESTRICT $WORKDIR/${PDBID}.btest $NTSTCNT $RFRRAT $RMSZB $RMSZA`

  #Print values
  echo " o Best B weight: $BBEST" | tee -a $LOG
else
  set BBEST = overall
  set BTEST = 0
endif

######################################################### re-refinement ##################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure re-refinement ******" | tee -a $LOG
echo "-Refining with geometric restraint weights: $WEIGHTS" | tee -a $LOG


#Run refmac with predefined matrix weights
foreach WGT (`echo $WEIGHTS`)

  #Return label for job launching
refirunning:

  #Only launch new jobs when the number of cores is not exceeded
  #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
  jobs > $WORKDIR/jobs.log

  if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_refin.pdb \
    XYZOUT $WORKDIR/${PDBID}_refmac${WGT}.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_refmac${WGT}.mtz \
    $LIBLIN \
    $TLSFILS \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${WGT}.log &
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      $BCMD
      refi type REST resi MLKF meth CGMAT bref $BREFTYPE
      $REFIRES
      tlsd waters exclude
      ncyc $NCYCLE
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG $WGTTYPE $WGT
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $OCCCMD
      $HARMCMD
      $METALCMD
      $LIBGCMD
      $RESTCMD
      temp $BBEST
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
      $ANOMCMD
      pdbout copy remarks 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      $KRABWGTCMD
      $KRABCMD
      END
eof
  else
    #Wait a bit to start again
    sleep 10
    goto refirunning
  endif
end

#Wait for the jobs to finish
wait

#Check for problems by seeing if the output mtz file exists
foreach WGT (`echo $WEIGHTS`)
  if (! -e $WORKDIR/${PDBID}_refmac$WGT.mtz) then
    echo " o Problem with refmac. Cannot continue." | tee -a $LOG
    echo "COMMENT: refmac: error in re-refinement" >> $WHYNOT
    echo "PDB-REDO,$PDBID"                         >> $WHYNOT
    if ($SERVER == 1) then
      #Write out status files
      touch $STDIR/stoppingProcess.txt
      touch $STDIR/processStopped.txt
    endif
    cd $BASE
    exit(1)
  endif
end

################################################# Find best results ######################################################

#Report
foreach WGT (`echo $WEIGHTS`)
  #Mine out the R-free
  set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_${WGT}.log | head -n 1 | awk '{print $3}'`
  echo " o Performed refinement with $WGTTYPE weight $WGT (R-free = $TFREE)" | tee -a $LOG
end


echo "-Selecting best geometric restraint weight" | tee -a $LOG

#No known errors in the refinement
set TLSERR  = 0

#Create second .refi file
if ($GOTR == 0 || $ZCALERR == 1) then
  echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}.refi
else
  echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}.refi
endif
foreach WGT (`echo $WEIGHTS`)
  if (-e $WORKDIR/${PDBID}_refmac${WGT}.pdb) then
    set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_${WGT}.log | head -n 1`
    echo "$WGT $LINE" >> $WORKDIR/${PDBID}.refi
  else
    echo "$WORKDIR/${PDBID}_refmac${WGT}.pdb is missing." | tee -a $LOG
    set TLSERR = 1
  endif
end

#Pick the best re-refined structure
set TLSBEST = `$TOOLS/picker $NEWMODEL $ESTRICT $WORKDIR/${PDBID}.refi $NTSTCNT $RFRRAT $RMSZB $RMSZA`

if ($TLSBEST == none) then

  #Set R(-free)
  set RTLS  = $RCAL
  set RFTLS = $RFCAL

  #Copy files
  cp $WORKDIR/${PDBID}_refin.pdb $WORKDIR/${PDBID}_besttls.pdb
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_ttest$OPTTLSG.mtz $WORKDIR/${PDBID}_besttls.mtz
    cp $WORKDIR/${PDBID}_ttest$OPTTLSG.log $WORKDIR/${PDBID}_besttls.log
  else
    cp $WORKDIR/${PDBID}_0cyc.mtz $WORKDIR/${PDBID}_besttls.mtz
    cp $WORKDIR/${PDBID}_0cyc.log $WORKDIR/${PDBID}_besttls.log
  endif

  #No use calculating values here
  set SIGRFTLS = 'NA'
  set RFTLSUNB = 'NA'
  set RFTLSZ   = 'NA'

else

  #Copy files
  cp $WORKDIR/${PDBID}_refmac${TLSBEST}.pdb $WORKDIR/${PDBID}_besttls.pdb
  cp $WORKDIR/${PDBID}_refmac${TLSBEST}.mtz $WORKDIR/${PDBID}_besttls.mtz
  cp $WORKDIR/${PDBID}_${TLSBEST}.log $WORKDIR/${PDBID}_besttls.log

  #Set R(-free)
  set RTLS  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | awk '{print $2}'`
  set RFTLS = `tail -n $LOGSTEP $WORKDIR/${PDBID}_besttls.log | head -n 1 | awk '{print $3}'`
  setenv PRTLS  $RTLS   #Set value for bits of Perl code in the script
  setenv PRFTLS $RFTLS  #Set value for bits of Perl code in the script

  #Validate R-free
  set SIGRFTLS = `perl -e 'printf ("%.4f\n", $ENV{PRFTLS}/sqrt($ENV{PNTSTCNT}));'`
  setenv PSIGRFTLS $SIGRFTLS
  set RFTLSUNB = `perl -e 'printf ("%.4f\n", $ENV{PRTLS}*$ENV{PRFRRAT});'`
  setenv PRFTLSUNB $RFTLSUNB
  set RFTLSZ   = `perl -e 'printf ("%.2f\n", ($ENV{PRFTLSUNB}-$ENV{PRFTLS})/$ENV{PSIGRFTLS});'`

endif
#Print values
echo " o Best geometric restraint weight: $TLSBEST" | tee -a $LOG

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Re-refinement details ******" | tee -a $LOG
echo "Best matrix weight: $TLSBEST" | tee -a $LOG
echo "Resulting R-factor: $RTLS"    | tee -a $LOG
echo "Resulting R-free  : $RFTLS"   | tee -a $LOG
echo "'Unbiased' R-free : $RFTLSUNB" | tee -a $LOG
echo "sigma(R-free)     : $SIGRFTLS" | tee -a $LOG
echo "R-free Z-score    : $RFTLSZ"   | tee -a $LOG

####################################################### Clean up round 1 #################################################

#Delete files...

#General
rm $WORKDIR/mtz_creation.log
rm $WORKDIR/raw.mtz
rm $WORKDIR/unique.mtz
rm $WORKDIR/merged.mtz
if (-e $WORKDIR/rawbu.mtz) then
  rm $WORKDIR/rawbu.mtz
endif
if (-e $WORKDIR/${PDBID}_0cycv1.log) then
    rm $WORKDIR/${PDBID}_0cycv1.log
    rm $WORKDIR/${PDBID}_0cycv1.pdb
endif

#TLS groups
if ($DOTLS == 1) then
  foreach TLSG (`ls ????.tls | cut -c 1-4`)
    rm $WORKDIR/${TLSG}_out.tls
    rm $WORKDIR/${PDBID}_ttest$TLSG.mtz
  end
endif

#B-weight
if ($BREFTYPE == ISOT || $BREFTYPE == ANISOT) then
  foreach BWGT (`echo $BWEIGHTS`)
    rm $WORKDIR/${PDBID}_btest$BWGT.pdb
    rm $WORKDIR/${PDBID}_btest$BWGT.mtz
    rm $WORKDIR/${PDBID}_btest$BWGT.log
  end
endif

#Re-refinement
foreach WGT (`echo $WEIGHTS`)
  rm $WORKDIR/${PDBID}_refmac${WGT}.pdb
  rm $WORKDIR/${PDBID}_refmac${WGT}.mtz
  rm $WORKDIR/${PDBID}_${WGT}.log
end

############################################   Validate structures  ######################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Validation details ******" | tee -a $LOG
echo "-Running WHAT_CHECK" | tee -a $LOG

#Initialise
set WCERR = 0

#Go to temporary running directory
setenv WCWORK $WORKDIR/wctemp
mkdir -p $WCWORK
cd $WCWORK

#Get the PDB file
cp $WORKDIR/${PDBID}_besttls.pdb $WCWORK

#Do the actual validation
$WC/bin/whatcheck $WCWORK/${PDBID}_besttls.pdb Y Y Y >& $WORKDIR/wc.log

#Check for an output file
if (-e $WCWORK/pdbout.txt) then
  #Do Nothing
else
  #Give warning
  echo " o WHAT_CHECK failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate re-refined structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                          >> $DEBUG
endif

#Create webpage
$WC/bin/pdbout2html >>& $WORKDIR/wc.log

#Check index.html completeness
if (`grep -c "Final summary" $WCWORK/pdbout.html` == 0 && `grep -c "Summary report" $WCWORK/pdbout.html` == 0) then
  echo " o WHAT_CHECK failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate re-refined structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                          >> $DEBUG
endif

#Clean up on aisle six
mkdir -p $WORKDIR/wc

if (-e $WCWORK/pdbout.txt) then
  mv $WCWORK/pdbout.txt $WORKDIR/wc/
endif
if (-e $WCWORK/check.db) then
  bzip2 $WCWORK/check.db
  mv $WCWORK/check.db.bz2 $WORKDIR/wc/
endif
#Server only
if ($SERVER == 1) then
  mv $WCWORK/pdbout.html $WORKDIR/wc/index.html
  mv $WCWORK/*.gif $WORKDIR/wc/ >& /dev/null
endif


#Go to start directory
cd $WORKDIR
rm -rf $WCWORK

#Extract stats from original structure
set ONATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/pdb${PDBID}.old`
if (-e $PDBOUT) then
  if (`grep -a -c '1st generation packing quality :' $PDBOUT` == 0) then
    set OZPAK1 = 'NA'
  else
    set OZPAK1 = `grep -a '1st generation packing quality :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $PDBOUT` == 0) then
    set OZPAK2 = 'NA'
  else
    set OZPAK2 = `grep -a '2nd generation packing quality :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $PDBOUT` == 0) then
    set OZRAMA = 'NA'
  else
    set OZRAMA = `grep -a 'Ramachandran plot appearance   :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $PDBOUT` == 0) then
    set OCHI12 = 'NA'
  else
    set OCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $PDBOUT` == 0) then
    set OBCONF = 'NA'
  else
    set OBCONF = `grep -a 'Backbone conformation          :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $PDBOUT` == 0) then
    set OBRMSZ = 'NA'
  else
    set OBRMSZ = `grep -a 'Bond lengths                   :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $PDBOUT` == 0) then
    set OARMSZ = 'NA'
  else
    set OARMSZ = `grep -a 'Bond angles                    :' $PDBOUT | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if ($OLDWC == 1) then
    set OWBMPS = 'NA'
    if (`grep -a -A1 -E '^.{31,32}<-->' $PDBOUT | tail -n 1 | grep -c 'total of'` == 0) then
      set OBUMPS = `grep -a -c -E '^.{31,32}<-->' $PDBOUT`
    else
      set OBUMPS = `grep -a -A1 -E '^.{31,32}<-->' $PDBOUT | tail -n 1 | awk '{print $8}'`
    endif
  else
    if (`grep -a -c 'Total number of bumps:' $PDBOUT` == 0) then
      set OBUMPS = 0
      set OSBMPL = 0
      set OWBMPS = 0.000
    else
      set OBUMPS = `grep -a 'Total number of bumps:' $PDBOUT | cut -c 24-28`
      set OSBMPL = `grep -a 'Total squared bump value:' $PDBOUT | cut -c 27-33`
      set OWBMPS = `echo $ONATOM $OSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
    endif
  endif
  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set OHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set OHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $PDBOUT | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set OHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set OHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $PDBOUT | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ OHBUNS = ($OHBDON + $OHBACC)
  if (`grep -a -c 'Buried donors:' $PDBOUT` == 0) then
    set OHBSAT = 'NA'
  else
    @ OBHBDA = (`grep -a 'Buried donors:' $PDBOUT | awk '{print $3}'` + `grep -a 'Buried acceptors:' $PDBOUT | awk '{print $3}'`)
    set OHBSA1 = `grep -a 'with a H-bond:' $PDBOUT | awk '{SUM += $5} END {print SUM}'`
    set OHBSA2 = `grep -a 'with a poor H-bond:' $PDBOUT | awk '{SUM += $6} END {print 0.5*SUM}'`
    set OHBSA3 = `grep -a 'with only a very poor H-bond:' $PDBOUT | awk '{SUM += $8} END {print 0.25*SUM}'`
    set OHBSA4 = `grep -a 'essentially without H-bond:' $PDBOUT | awk '{SUM += $5} END {print 0.125*SUM}'`
    set OHBSAT = `echo $OHBSA1 $OHBSA2 $OHBSA3 $OHBSA4 $OBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
  endif
else
  set OZPAK1 = 'NA'
  set OZPAK2 = 'NA'
  set OZRAMA = 'NA'
  set OCHI12 = 'NA'
  set OBCONF = 'NA'
  set OBRMSZ = 'NA'
  set OARMSZ = 'NA'
  set OBUMPS = 'NA'
  set OWBMPS = 'NA'
  set OHBUNS = 'NA'
  set OHBSAT = 'NA'
endif

#Extract stats from re-refined structure
set NNATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_besttls.pdb`
if (-e $WORKDIR/wc/pdbout.txt) then
  if (`grep -a -c '1st generation packing quality :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZPAK1 = 'NA'
  else
    set NZPAK1 = `grep -a '1st generation packing quality :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZPAK2 = 'NA'
  else
    set NZPAK2 = `grep -a '2nd generation packing quality :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NZRAMA = 'NA'
  else
    set NZRAMA = `grep -a 'Ramachandran plot appearance   :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NCHI12 = 'NA'
  else
    set NCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NBCONF = 'NA'
  else
    set NBCONF = `grep -a 'Backbone conformation          :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NBRMSZ = 'NA'
  else
    set NBRMSZ = `grep -a 'Bond lengths                   :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $WORKDIR/wc/pdbout.txt` == 0) then
    set NARMSZ = 'NA'
  else
    set NARMSZ = `grep -a 'Bond angles                    :' $WORKDIR/wc/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if ($OLDWC == 1) then
    set NWBMPS = 'NA'
    if (`grep -a -A1 -E '^.{31,32}<-->' $WORKDIR/wc/pdbout.txt | tail -n 1 | grep -c 'total of'` == 0) then
      set NBUMPS = `grep -a -c -E '^.{31,32}<-->' $WORKDIR/wc/pdbout.txt`
    else
      set NBUMPS = `grep -a -A1 -E '^.{31,32}<-->' $WORKDIR/wc/pdbout.txt | tail -n 1 | awk '{print $8}'`
    endif
  else
    if (`grep -a -c 'Total number of bumps:' $WORKDIR/wc/pdbout.txt` == 0) then
      set NBUMPS = 0
      set NSBMPL = 0
      set NWBMPS = 0.000
    else
      set NBUMPS = `grep -a 'Total number of bumps:' $WORKDIR/wc/pdbout.txt | cut -c 24-28`
      set NSBMPL = `grep -a 'Total squared bump value:' $WORKDIR/wc/pdbout.txt | cut -c 27-33`
      set NWBMPS = `echo $NNATOM $NSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
    endif
  endif
  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set NHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set NHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set NHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set NHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wc/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ NHBUNS = ($NHBDON + $NHBACC)
  if (`grep -a -c 'Buried donors:' $WORKDIR/wc/pdbout.txt` == 0) then
    set NHBSAT = 'NA'
  else
    @ NBHBDA = (`grep -a 'Buried donors:' $WORKDIR/wc/pdbout.txt | awk '{print $3}'` + `grep -a 'Buried acceptors:' $WORKDIR/wc/pdbout.txt | awk '{print $3}'`)
    set NHBSA1 = `grep -a 'with a H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $5} END {print SUM}'`
    set NHBSA2 = `grep -a 'with a poor H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $6} END {print 0.5*SUM}'`
    set NHBSA3 = `grep -a 'with only a very poor H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $8} END {print 0.25*SUM}'`
    set NHBSA4 = `grep -a 'essentially without H-bond:' $WORKDIR/wc/pdbout.txt | awk '{SUM += $5} END {print 0.125*SUM}'`
    set NHBSAT = `echo $NHBSA1 $NHBSA2 $NHBSA3 $NHBSA4 $NBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
  endif
else
  set NZPAK1 = 'NA'
  set NZPAK2 = 'NA'
  set NZRAMA = 'NA'
  set NCHI12 = 'NA'
  set NBCONF = 'NA'
  set NBRMSZ = 'NA'
  set NARMSZ = 'NA'
  set NBUMPS = 'NA'
  set NWBMPS = 'NA'
  set NHBUNS = 'NA'
  set NHBSAT = 'NA'
endif

#Run FoldX (if available)
if ($?FOLDX) then
  echo "-Running FoldX" | tee -a $LOG
  #Copy over the rotabase.txt file
  cp $ROTABASE $WORKDIR/rotabase.txt

  #Run the original PDB file
  $FOLDX --command=Stability --pdb-dir=$WORKDIR --pdb=${PDBID}_0cyc.pdb --output-file=$PDBID >>& $WORKDIR/foldx.log

  #Set the value
  if (`grep -c 'Total          =' $WORKDIR/foldx.log` != 0) then
    set OGFOLD = `grep 'Total          =' $WORKDIR/foldx.log | tail -n 1 | awk '{print $3}'`
  else
    set OGFOLD = 'NA'
  endif

  #Run the rerefined PDB file
  $FOLDX --command=Stability  --pdb-dir=$WORKDIR --pdb=${PDBID}_besttls.pdb --output-file=$PDBID >>& $WORKDIR/foldx.log

  #Set the value
  if (`grep -c 'Total          =' $WORKDIR/foldx.log` != 0) then
    set NGFOLD = `grep 'Total          =' $WORKDIR/foldx.log | tail -n 1 | awk '{print $3}'`
  else
    set NGFOLD = 'NA'
  endif

else
  set OGFOLD = 'NA'
  set NGFOLD = 'NA'
endif

#Run distel if needed
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  #Consolidate the distance restraints
  if (-e $WORKDIR/homology.rest) then
    cat $WORKDIR/homology.rest > $WORKDIR/allhb.rest
  endif
  if (-e $WORKDIR/hbond.rest) then
    cat $WORKDIR/hbond.rest >> $WORKDIR/allhb.rest
  endif
  
  #Calculate rmsZ values. PROGRAM: distel.py
  set OHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set NHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`


  #Clean up
  rm $WORKDIR/allhb.rest
else 
  set OHRMSZ = 'NA'
  set NHRMSZ = 'NA'
endif

#Run distel for torsion angle restraints
if ($NKRABREST > 0) then
  set OKRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/krab.rest | tail -n 1 | awk '{print $5}'`
  set NKRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/krab.rest | tail -n 1 | awk '{print $5}'`
else
  set OKRMSZ = 'NA'
  set NKRMSZ = 'NA'
endif

#Print results
echo '                                     Before  After' | tee -a $LOG
echo "1st generation packing quality     : $OZPAK1 $NZPAK1" | tee -a $LOG
echo "2nd generation packing quality     : $OZPAK2 $NZPAK2" | tee -a $LOG
echo "Ramachandran plot appearance       : $OZRAMA $NZRAMA" | tee -a $LOG
echo "chi-1/chi-2 rotamer normality      : $OCHI12 $NCHI12" | tee -a $LOG
echo "Backbone conformation              : $OBCONF $NBCONF" | tee -a $LOG
echo " " | tee -a $LOG
echo "Bond length RMS Z-score            : $OBRMSZ $NBRMSZ" | tee -a $LOG
echo "Bond angle RMS Z-score             : $OARMSZ $NARMSZ" | tee -a $LOG
echo " " | tee -a $LOG
echo "Total number of bumps              : $OBUMPS $NBUMPS" | tee -a $LOG
echo "Weighted bump severity score       : $OWBMPS $NWBMPS" | tee -a $LOG
echo " " | tee -a $LOG
echo "Unsatisfied H-bond donors/acceptors: $OHBUNS $NHBUNS" | tee -a $LOG
echo "H-bond satisfaction fraction       : $OHBSAT $NHBSAT" | tee -a $LOG
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  echo "H-bond restraint RMS Z-score       : $OHRMSZ $NHRMSZ" | tee -a $LOG
endif
if ($NKRABREST > 0) then
  echo "KRAB restraint RMS Z-score         : $OKRMSZ $NKRMSZ" | tee -a $LOG
endif
if ($?FOLDX) then
  echo " " | tee -a $LOG
  echo "Gibbs folding energy (kcal/mol)    : $OGFOLD $NGFOLD" | tee -a $LOG
endif

################################################## Rebuild structure #####################################################

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Structure rebuilding ******" | tee -a $LOG

#Temporary HETATM and LINK workaround
mv $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_besttls.old
sed '/MSE/s/HETATM/ATOM  /g' $WORKDIR/${PDBID}_besttls.old > $WORKDIR/${PDBID}_besttls.mse
$TOOLS/stripper -v $SMODE \
$WORKDIR/${PDBID}_besttls.mse \
$WORKDIR/${PDBID}_besttls.pdb \
$TOOLS/pdb_redo.dat \
>> $WORKDIR/stripper.log

#Set fall-back values
set NWATDEL = 0
set NBBFLIP = 0
set NSCBLT  = 0
set NSCFLIP = 0
set NCHIRFX = 0
set BUILT   = 0  #The model was not rebuilt
touch $WORKDIR/${PDBID}_pepflip.log #We need a pepflip log file for coot_tour

if ($DOREBUILD == 0) then
  echo "-All rebuilding steps are skipped" | tee -a $LOG
  cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_built.pdb
else
  #Calculate BBUILD
  set BBUILD = `perl -e 'printf ("%.2f\n", 7*log($ENV{BSET}));'`

  #Run dssp to find secondary structure elements
  if ($GOT_PROT == T) then
    echo "-Assigning secondary structure with DSSP" | tee -a $LOG
    mkdssp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/$PDBID.dssp >& $WORKDIR/dssp.log
    if($status) then
      echo " o DSSP failed" | tee -a $LOG
      echo "   * Not using secondary structure in rebuilding" | tee -a $LOG
      echo "COMMENT: DSSP: general error (1)" >> $DEBUG
      echo "PDB-REDO,$PDBID"                  >> $DEBUG
    endif

    #Is there usable output?
    if (-e $WORKDIR/$PDBID.dssp) then
      set DSSPFILE = $WORKDIR/$PDBID.dssp
    else
      set DSSPFILE = ""
    endif
  endif

  #Check the model completeness and create an omit map of sorts if needed
  if (`perl -e 'if ($ENV{COMPLETED} < 50.0) {print "1"} else {print "0"}'` == 1) then
    #Make a map with SER/PRO/VAL/THR side chains trimmed
    grep -v '^A[NT][IO][SM].........[OC][GD][ 12].[SVTP][EAHR][RLO]' $WORKDIR/${PDBID}_besttls.pdb > $WORKDIR/${PDBID}_besttls_trim.pdb

    #Run Refmac 0cyc to make the maps
    refmac5 \
    XYZIN  $WORKDIR/${PDBID}_besttls_trim.pdb \
    XYZOUT $WORKDIR/${PDBID}_scomit.pdb \
    HKLIN  $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/${PDBID}_scomit.mtz \
    $LIBLIN \
    $TLSFILS \
    $SCATLIN \
<<eof >$WORKDIR/${PDBID}_${WGT}.log
      $SCATTERCMD
      make check NONE
      make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
        ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
      refi type REST resi MLKF meth CGMAT bref $BREFTYPE
      $REFIRES
      tlsd waters exclude
      ncyc 0
      scal type $SOLVENT $SCALING
      solvent YES
      $MASKPAR
      $LOWMEM
      weight $WGTSIG $WGTTYPE $WGT
      monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
        chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
      $TWIN
      $NCSTYPE
      $NCSALIGN
      $NCSNEIGH
      $NCSSTRICT
      $JELLY
      $TORSION
      $OCCCMD
      $HARMCMD
      $METALCMD
      $LIBGCMD
      $RESTCMD
      temp $BBEST
      blim 2.0 999.0
      labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
      pdbout copy remarks 280 350
      NOHARVEST
      $HOMOLWGTCMD
      $HOMOLCMD
      $HBONDWGTCMD
      $HBONDCMD
      $KRABWGTCMD
      $KRABCMD
      END
eof

#     cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_besttls_noomit.mtz
#     comit \
#     --mtzin  $WORKDIR/${PDBID}_besttls_noomit.mtz \
#     --mtzout $WORKDIR/${PDBID}_besttls.mtz \
#     --colin-fo FP,SIGFP \
#     --colin-fc FC_ALL,PHIC_ALL > $WORKDIR/comit.log
#
#     #Set map coefficients for building
#     set BUILDF = 'omit.F_phi.F'
#     set BUILDP = 'omit.F_phi.phi'
    #Set map coefficients for building
    set SCBUILDMTZ = $WORKDIR/${PDBID}_scomit.mtz
    set BUILDF     = 'FWT'
    set BUILDP     = 'PHWT'
  else
    #Set map coefficients for building
    set SCBUILDMTZ = $WORKDIR/${PDBID}_besttls.mtz
    set BUILDF     = 'FWT'
    set BUILDP     = 'PHWT'

  endif


  #Count the number of waters
  set NWATER = `grep -c -E '^[AH][TE][OT][MA].{13}HOH' $WORKDIR/${PDBID}_besttls.pdb`

  if ($DOCENTRIFUGE == 1 && $NWATER > 0) then
    #Delete poor waters
    echo "-Removing waters without density" | tee -a $LOG

    $TOOLS/centrifuge \
    <<eof >& $WORKDIR/${PDBID}_centrifuge.log
      #Settings
      InputPDB = "$WORKDIR/${PDBID}_besttls.pdb"
      PDBOutputFilename = "$WORKDIR/${PDBID}_centrifuge.pdb"
      MessageFilename = "$WORKDIR/${PDBID}_centrifuge.msg"
      XMLOutputFilename = "$WORKDIR/${PDBID}_centrifuge.xml"

      #residue names that represent waters. Default HOH WAT H2O EAU
      WaterNames = HOH
      #Recognised waters, and dummies, with a density fit lower than this threshold will be removed
      WaterRejectionThreshold = 0.30
      #If true, everything except Dummies and recognised Waters is placed in the density map, default False
      PlaceNonSolventInDensity = False
      #Colon separated list of water ids which won't be affected by the methods of centrifuge
      UnaffectedWaterList = "$H2O_KEEP"
      #InputMap or here: input mtz
      InputMTZ = "$WORKDIR/${PDBID}_besttls.mtz"
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      #type of method used to determine the density fit
      CFitTarget::Type = accelerated
      UniformAtomBFactor = "$BBUILD"
      UniformAtomRadius = 0.74
      UseAtomicB = True

      #Factor with which to invert the density of placed atoms
      AntiBumpFactor = 2.0

      #Program information
      ProgramName = centrifuge
      MessageLevel = 6
      AbortLevel = 8

      KeepPDBheaderInfo = true
      StoreOriginalChainAndSegID = true
      SelectByOrigChainID = true
      KeepSideChain = true
      KeepFragmentLongerThan = 0
      TrustWaters = false

      RemoveBasedOnDensity = true
eof
      if($status) then
        #Write general debug statement
        echo " o Problem with centrifuge" | tee -a $LOG
        echo "COMMENT: error in centrifuge" >> $DEBUG
        echo "PDB-REDO,$PDBID"              >> $DEBUG
        cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_centrifuge.pdb
      else
        #Write debug messsages for building tools
        if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_centrifuge.msg` > 0) then
          grep -H WRNG $WORKDIR/${PDBID}_centrifuge.msg >> $DEBUGB
        endif
      endif

    set NWATDEL = `grep -c 'removing atom' $WORKDIR/${PDBID}_centrifuge.log`

    #If waters were deleted, consider model rebuilt
    if ($NWATDEL > 0) then
      set BUILT = 1
      #Delete external restraints refering to the waters
      if ($RESTCMD != "") then
        foreach WAT (`grep 'removing atom' $WORKDIR/${PDBID}_centrifuge.log | cut -c 16-20`)
          cp $WORKDIR/external.rest $WORKDIR/external.bak
          set CHID = `echo $WAT | cut -c 1-1`
          set RESN = `echo $WAT | cut -c 2-`
          grep -v "chain $CHID resi $RESN" $WORKDIR/external.bak > $WORKDIR/external.rest
        end
      endif
    endif

  else
    echo "-Skipping centrifuge run" | tee -a $LOG
    cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_centrifuge.pdb
  endif

  #Further rebuilding steps

  #Perform peptide flips if resolution is adequate and there is protein
  if ($DOPEPFLIP == 1 && $GOT_PROT == T && `perl -e 'if ($ENV{URESO} < 3.30) {print "1"}'` == 1) then
    echo "-Performing peptide flips" | tee -a $LOG

    #Run pepflip
    $TOOLS/pepflip \
<<eof >& $WORKDIR/${PDBID}_pepflip.log
      ProgramName = pepflip
      MessageFilename = $WORKDIR/${PDBID}_pepflip.msg
      MessageLevel = 6
      XMLOutputFilename = $WORKDIR/${PDBID}_pepflip.xml
      AbortLevel = 8

      #input pdb, and its handling
      InputPDB = $WORKDIR/${PDBID}_centrifuge.pdb
      SelectByOrigChainID =true
      KeepPDBheaderInfo =true

      #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
      InputMTZ = $WORKDIR/${PDBID}_besttls.mtz
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      #Use this masked map
      OutputMaskedNonAAmap = $WORKDIR/masked.map

      #Output map with the pepflip version of asu limits, needed for loopfit
      #OutputExtendedMap = $WORKDIR/${PDBID}_map.ext
      #density settings for the computations on the difference density at the oxygen
      DifDensity::InputMTZ = $WORKDIR/${PDBID}_besttls.mtz
      DifDensity::FWTLabel = DELFWT
      DifDensity::PHIWTLabel = PHDELWT
      DifDensity::CFitTarget::Type = accelerated
      DifDensity::CFitTarget::ConvertToZscore = false

      UniformAtomRadius= 0.74
      AntiBumpFactor= 2.0
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"
      UniformAtomBFactor = "$BBUILD"

      FlipMode = repair
      #Output pdb with the flipped peptides
      FlippedOutputPDB = $WORKDIR/${PDBID}_pepflip.pdb

      #Exclude lists
      DontFlipListN = "$BBN_KEEP"
      DontFlipListO = "$BBO_KEEP"

      #Optional filename (output from DSSP) to set the secondary structure of the complex
      DSSPfilename = $DSSPFILE
      #Number of residues (default 2) at the ends of secondary structure, which won't be trusted
      BufferSS=1
      #If true (default false), and the secondary structure is set, trust the middle of helices (H)
      #and don't consider those residues for flips
      TrustHelices = true
      Trust310Helices = true
      TrustPiHelices = true
      TrustBetaStrands = true

      #Refinement parameters
      Loopfit::UseDefaultLoopfit =false
      Loopfit::UseMiniRSR=true
      Loopfit::InputMap= $WORKDIR/masked.map
      Loopfit::LoopfitExe = coot-mini-rsr
      Loopfit::LoopfitLog = $WORKDIR/pepflip.fit
      Loopfit::MiniRSRtorsions = true
      Loopfit::MiniRSRrama = true
      Loopfit::MiniRSRweight = 10.00
      #temporary pdb file for loopfit/mini-rsr
      Loopfit::TmpOutputPDB = $WORKDIR/${PDBID}_fit.pdb
      #Location of the CCP4 dictionary
      Loopfit::MonomerLib = $CLIBD_MON
      WeightGeoGooF = 24

      #local dictionary files
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      #Methods for density computations
      CFitTarget::Type = cubic
      CFitTarget::MaskRadius = 2.0
      CFitTarget::UseObservedDensityAroundAtoms = "CA CB"

      #In repair mode - If true, place the part of the flanking residues that doesn't belong to the peptide in the map,
      #default false
      RepairPlaceFlanking = true
      #If true, place known residues in the map, before checking the density fit of the peptide - This will slow down
      #the progress significantly, default false
      CheckDensityPlaceKnown = false

      #Reject flip if the new density < average(density) + 'MinSigmaLevel'*sigma(density)
      MinSigmaLevel = -3.5
      #Ratio between the density fit of the original peptide and the flipped version, to be considered as possible flip (default 0.9)
      MinRatioOrigFlip =0.9
      #Density fit of Oxygen in the difference map, indicating a clear positive monopole (default 2.)
      DensityPosOmonopole = 2.4
      #Density fit of Oxygen in the difference map, indicating a clear negative monopole (default -3.)
      DensityNegOmonopole = -2.4
eof
    if($status) then
      echo " o Problem with pepflip. Cannot continue." | tee -a $LOG
      echo "COMMENT: error in pepflip" >> $WHYNOT
      echo "PDB-REDO,$PDBID"           >> $WHYNOT
      if ($SERVER == 1) then
        #Write out status files
        touch $STDIR/stoppingProcess.txt
        touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      #Write debug messsages for building tools
      if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_pepflip.msg` > 0) then
        grep -H WRNG $WORKDIR/${PDBID}_pepflip.msg >> $DEBUGB
      endif
    endif

    #Count the number of flips
    set NBBFLIP = `grep -A 80 'List of peptides' $WORKDIR/${PDBID}_pepflip.log | grep -c -E '^.[0-9]'`

    #If peptides were flipped, consider model rebuilt
    if ($NBBFLIP > 0) then
      set BUILT = 1
    endif

  else
    echo "-Skipping pepflip run"
    cp $WORKDIR/${PDBID}_centrifuge.pdb $WORKDIR/${PDBID}_pepflip.pdb
  endif

  #Rebuild side chains if resolution is adequate and there is protein
  if ($DOSCBUILD == 1 && $GOT_PROT == T && `perl -e 'if ($ENV{URESO} < 3.30) {print "1"}'` == 1) then
    echo "-Refitting side chains" | tee -a $LOG

    #Run dssp again to find secondary structure elements if there were any peptide flips
    if ($NBBFLIP > 0) then
      mkdssp $WORKDIR/${PDBID}_centrifuge.pdb $WORKDIR/$PDBID.dssp >>& $WORKDIR/dssp.log
      if($status) then
        echo " o Problem with DSSP" | tee -a $LOG
        echo "   * Not using secondary structure in SideAide" | tee -a $LOG
        echo "COMMENT: DSSP: general error (2)" >> $DEBUG
        echo "PDB-REDO,$PDBID"                  >> $DEBUG
      endif

      #Is there usable output?
      if (-e $WORKDIR/$PDBID.dssp) then
        set DSSPFILE = $WORKDIR/$PDBID.dssp
      else
        set DSSPFILE = ""
      endif
    endif

    #Run SideAide
    $TOOLS/SideAide \
  <<eof >& $WORKDIR/${PDBID}_scbuild.log
      ProgramName = SideAide
      MessageFilename = $WORKDIR/${PDBID}_scbuild.msg
      MessageLevel = 6
      XMLOutputFilename = $WORKDIR/${PDBID}_scbuild.xml
      AbortLevel = 8

      #Input pdb, and its handling
      #IncludeChains = '-'
      InputPDB = "$WORKDIR/${PDBID}_pepflip.pdb"

      #list of residue names representing waters, default HOH H2O WAT EAU
      #WaterNames = HOH WAT H2O EAU

      #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
      InputMTZ = "$SCBUILDMTZ"
      FWTLabel = $BUILDF
      PHIWTLabel = $BUILDP
      SpaceGroup = "$SPACEGROUP"
      XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

      #B factor used in map calculation (for every atom in every residue):
      UniformAtomBFactor = $BBUILD
      UniformAtomRadius = 0.74
      AntiBumpFactor = 2.0

      #Optional filename (output from DSSP) to set the secondary structure of the complex
      DSSPfilename = $DSSPFILE

      #Output pdb
      PDBOutputFilename = "$WORKDIR/${PDBID}_scbuild.pdb"

      #selection of the side chains to rebuild/flip
      #list of single residues:
      RebuildAll = true #If false, no residues are rebuilt (use for flipping without rebuilding)
      UseAvBforNewSCatoms = true #default
      #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false

      #RebuildList = ":A33 :A66 :A69 :A87 :A167 :A245 :A275 :A304 :A314 :A412 :"
      DontRebuildList = "$NO_REBUILD"
      TrustedWaterList = "$H2O_KEEP"
      #list of residue regions
      #RebuildDefinition = "A36(7)A42:AAAA90(6)AAAA95"
      #list of residues to flip
      #FlipResidueList = ":A23 :"

      #local dictionary files
      DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
      SymmetryFilename = "$CLIBD/syminfo.lib"

      # Methods for density computations
      CFitTarget::Type = accelerated
      CFitTarget::MaskRadius = 2.0
      CFitTarget::UseObservedDensityAroundAtoms = "CA"
      RefineScope::CFitTarget::Type = accelerated
      FitRotamerScope::CFitTarget::Type = correlation

      #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
      DummyRejectionThreshold = -0.1
      RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

      #Advanced parameter settings
      KeepPDBheaderInfo = true
      KeepUnknownAsIs = true

      #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
      # rotamer fit exceeds a threshold (these will be fit and refined in a second round
      OrderOnRotamerValidation = true
      #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
      ValidationThresholds = "all 0.4 ALA 0.3"

      #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
      #Only used for CDummyResidues and CFullRes, i.e. ignored residues
      StoreOriginalChainAndSegID = true
      KeepSideChain = true
      KeepFragmentLongerThan = 0
      #TrustWaters = false
      AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MSE MET LYS ARG"
      #place residues outside IncludeChains in the map:
      InitTargetWithIgnored = true
      ThresholdNewRefinedSC = -0.04 #Only keep sidechains if the fit improves by this value

      #Parameters to shift the CA to try and find better rotamers, optional
      #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
      #Radius of the CA shell of shifted options
      ShiftingCAshell = 0.2
      #number of CA, or closest fibonaccinumber, in the shell of CA options
      ShiftingCAnumber = 8
      #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
      ShiftingCAthreshold = 10.
      #only keep a rotamer with a shifted CA if the density score is higher than this value
      ShiftingCAkeepThreshold = -10
      #-only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
      ShiftingCAimprovementThreshold = 0.05

      #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
      SelectByOrigChainID = true
      Overlap = 0
eof
    if($status) then
      echo " o Problem with SideAide. Cannot continue." | tee -a $LOG
      echo "COMMENT: error in SideAide" >> $WHYNOT
      echo "PDB-REDO,$PDBID"            >> $WHYNOT
      if ($SERVER == 1) then
	#Write out status files
	touch $STDIR/stoppingProcess.txt
	touch $STDIR/processStopped.txt
      endif
      cd $BASE
      exit(1)
    else
      #Write debug messsages for building tools
      if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_scbuild.msg` > 0) then
        grep -H WRNG $WORKDIR/${PDBID}_scbuild.msg >> $DEBUGB
      endif
    endif

    #Count the number of completed side chains.
    set NSCBLT  = `grep -c 'building the complete side chain' $WORKDIR/${PDBID}_scbuild.log`

    #Some side chains must have changed, so consider the model rebuilt
    set BUILT = 1

  else
    echo "-Skipping side chain rebuilding" | tee -a $LOG
    cp $WORKDIR/${PDBID}_pepflip.pdb $WORKDIR/${PDBID}_scbuild.pdb
  endif

######################################## Flip side chains to optimise hydrogen bonding ###################################

  #Always do this if there is protein
  if ($GOT_PROT == T) then

    echo "-Validation-based rebuilding" | tee -a $LOG

    #Check for candidates using WHAT_CHECK
    #Go to temporary running directory
    mkdir -p $WORKDIR/flip
    cd $WORKDIR/flip

    #Get the PDB file
    cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/flip

    #Do the actual validation
    $WC/bin/whatcheck $WORKDIR/${PDBID}_scbuild.pdb Y Y Y >& $WORKDIR/flip/flip.log

    #Check for an output file
    if (-e $WORKDIR/flip/pdbout.txt) then

      #Create a fliplist
      $TOOLS/what_todo $WORKDIR/flip/pdbout.txt $WORKDIR/$PDBID.todo $WORKDIR/$PDBID.extracted

      #Do the flipping
      cd $WORKDIR
      set FLIPLIST = `head -n 2  $WORKDIR/$PDBID.todo | tail -n 1`
      echo " o Flipping side chains"

      #Run SideAide
      $TOOLS/SideAide \
<<eof >& $WORKDIR/${PDBID}_scflip.log
        ProgramName = SideAide
        MessageFilename = $WORKDIR/${PDBID}_scflip.msg
        MessageLevel = 6
        XMLOutputFilename = $WORKDIR/${PDBID}_scflip.xml
        AbortLevel = 8

        #Input pdb, and its handling
        InputPDB = "$WORKDIR/${PDBID}_scbuild.pdb"

        #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
        InputMTZ = "$WORKDIR/${PDBID}_besttls.mtz"
        FWTLabel = $BUILDF
        PHIWTLabel = $BUILDP
        SpaceGroup = "$SPACEGROUP"
        XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

        #B factor used in map calculation (for every atom in every residue):
        UniformAtomBFactor = $BBUILD
        UniformAtomRadius = 0.74
        AntiBumpFactor = 2.0

        #Optional filename (output from DSSP) to set the secondary structure of the complex
        DSSPfilename = $DSSPFILE

        #Output pdb
        PDBOutputFilename = "$WORKDIR/${PDBID}_built.pdb"

        #selection of the side chains to rebuild/flip
        #list of single residues:
        RebuildAll = false #If false, no residues are rebuilt (use for flipping without rebuilding)
        UseAvBforNewSCatoms = true #default
        #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false
        DontRebuildList = "$NO_REBUILD"
        TrustedWaterList = "$H2O_KEEP"
        #list of residues to flip
        FlipResidueList = "$FLIPLIST"

        #local dictionary files
        DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
        SymmetryFilename = "$CLIBD/syminfo.lib"

        #Methods for density computations
        CFitTarget::Type = accelerated
        CFitTarget::MaskRadius = 2.0
        CFitTarget::UseObservedDensityAroundAtoms = "CA"
        RefineScope::CFitTarget::Type = accelerated
        FitRotamerScope::CFitTarget::Type = correlation

        #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
        DummyRejectionThreshold = -0.1
        RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

        #Advanced parameter settings
        KeepPDBheaderInfo = true
        KeepUnknownAsIs = true

        #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
        # rotamer fit exceeds a threshold (these will be fit and refined in a second round
        OrderOnRotamerValidation = true
        #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
        ValidationThresholds = "all 0.4 ALA 0.3"

        #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
        #Only used for CDummyResidues and CFullRes, i.e. ignored residues
        StoreOriginalChainAndSegID = true
        KeepSideChain = true
        KeepFragmentLongerThan = 0
        #TrustWaters = false
        AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MET LYS ARG"
        #place residues outside IncludeChains in the map:
        InitTargetWithIgnored = true
        ThresholdNewRefinedSC = -0.04 #Only keep sidechains if the fit improves by this value

        #Parameters to shift the CA to try and find better rotamers, optional
        #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
        #Radius of the CA shell of shifted options
        ShiftingCAshell = 0.2
        #number of CA, or closest fibonaccinumber, in the shell of CA options
        ShiftingCAnumber = 8
        #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
        ShiftingCAthreshold = 10.
        #only keep a rotamer with a shifted CA if the density score is higher than this value
        ShiftingCAkeepThreshold = -10
        #only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
        ShiftingCAimprovementThreshold = 0.10

        #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
        SelectByOrigChainID = true
        Overlap = 0
eof
      if($status) then
        echo "   * Problem with SideAide. Cannot continue." | tee -a $LOG
        echo "COMMENT: error in side chain flipping" >> $WHYNOT
        echo "PDB-REDO,$PDBID"                       >> $WHYNOT
        if ($SERVER == 1) then
          #Write out status files
          touch $STDIR/stoppingProcess.txt
          touch $STDIR/processStopped.txt
        endif
        cd $BASE
        exit(1)
      else
        #Write debug messsages for building tools
        if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_scflip.msg` > 0) then
          grep -H WRNG $WORKDIR/${PDBID}_scflip.msg >> $DEBUGB
        endif
      endif

      #Get the number of flipped residues
      set NSCFLIP = `grep -c 'Flipping residue' $WORKDIR/${PDBID}_scflip.log`

      #If sidechains were flipped, consider model rebuilt
      if ($NSCFLIP > 0) then
        set BUILT = 1
      endif

      #Force rebuilding seriously problematic side chains
      if (`grep -c '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo` != 0) then
        set CHIRLIST = `grep -A 1 '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo | tail -n 1`
        cp $WORKDIR/${PDBID}_built.pdb $WORKDIR/${PDBID}_chirfix_in.pdb

        #Add a debug warning
        echo " o Chirality fixes are needed" | tee -a $LOG

        #Run SideAide
        $TOOLS/SideAide \
<<eof >& $WORKDIR/${PDBID}_chirfix.log
        ProgramName = SideAide
        MessageFilename = $WORKDIR/${PDBID}_chirfix.msg
        MessageLevel = 6
        XMLOutputFilename = $WORKDIR/${PDBID}_chirfix.xml
        AbortLevel = 8

        #Input pdb, and its handling
        #IncludeChains = '-'
        InputPDB = "$WORKDIR/${PDBID}_chirfix_in.pdb"

        #list of residue names representing waters, default HOH H2O WAT EAU
        #WaterNames = HOH WAT H2O EAU

        #Density information: Use either InputMap, or InputMTZ plus FWTLabel and PHIWTLabel
        InputMTZ = "$WORKDIR/${PDBID}_besttls.mtz"
        FWTLabel = $BUILDF
        PHIWTLabel = $BUILDP
        SpaceGroup = "$SPACEGROUP"
        XtalCell = "$AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA"

        #B factor used in map calculation (for every atom in every residue):
        UniformAtomBFactor = $BBUILD
        UniformAtomRadius = 0.74
        AntiBumpFactor = 2.0

        #Optional filename (output from DSSP) to set the secondary structure of the complex
        DSSPfilename = $DSSPFILE

        #Output pdb
        PDBOutputFilename = "$WORKDIR/${PDBID}_built.pdb"

        #selection of the side chains to rebuild/flip
        #list of single residues:
        RebuildAll = false #If false, no residues are rebuilt (use for flipping without rebuilding)
        UseAvBforNewSCatoms = true #default
        #BFactorForNewSCatoms = 15 #ouse if UseAvBforNewSCatoms = false

        RebuildList = "$CHIRLIST"
        TrustedWaterList = "$H2O_KEEP"
        #list of residue regions
        #list of residues to flip

        #local dictionary files
        DictionaryFilename = "$TOOLS/AArotaSS.XYZ"
        SymmetryFilename = "$CLIBD/syminfo.lib"

        #Methods for density computations
        CFitTarget::Type = accelerated
        CFitTarget::MaskRadius = 2.0
        CFitTarget::UseObservedDensityAroundAtoms = "CA"
        RefineScope::CFitTarget::Type = accelerated
        FitRotamerScope::CFitTarget::Type = correlation

        #Dummies and dubious ignored residues, i.e. possible waters, are removed when their density lies below this value
        DummyRejectionThreshold = -0.1
        RejectionScope::CFitTarget::type = accelerated #kies de fittarget die hoort by DummyRejectionThreshold.

        #if true, the order of fitting rotamers and refinement is not only based on the AA type, but only on whether the
        # rotamer fit exceeds a threshold (these will be fit and refined in a second round
        OrderOnRotamerValidation = true
        #list of AA type and their validation threshold ('all' can be used to set the threshold of every aa type).
        ValidationThresholds = "all 0.4 ALA 0.3"

        #Advanced parameter settings
        KeepPDBheaderInfo = true
        KeepUnknownAsIs = true

        #Below:  default is false, if true, the original chain ID and Seg ID are saved in the output pdb.
        #Only used for CDummyResidues and CFullRes, i.e. ignored residues
        StoreOriginalChainAndSegID = true
        KeepSideChain = true
        KeepFragmentLongerThan = 0
        #TrustWaters = false
        AminoAcidFittingOrder = "GLY ALA PRO VAL THR SER CYS TRP PHE TYR HIS ILE LEU ASP ASN GLU GLN MET LYS ARG"
        #place residues outside IncludeChains in the map:
        InitTargetWithIgnored = true
        ThresholdNewRefinedSC = -1.00 #Very low value to ensure that original side chains are not kept

        #Parameters to shift the CA to try and find better rotamers, optional
        #Values to twist the N-CA-CB angle over, to try and achieve better rotamer fitting;
        #Radius of the CA shell of shifted options
        ShiftingCAshell = 0.2
        #number of CA, or closest fibonaccinumber, in the shell of CA options
        ShiftingCAnumber = 8
        #only try to find a better rotamer if the density score is lower than this value, 10 = NO THRESHOLD
        ShiftingCAthreshold = 10.
        #only keep a rotamer with a shifted CA if the density score is higher than this value
        ShiftingCAkeepThreshold = -10
        #-only keep a rotamer with a shifted CA if the relative improvement in the score is higher than this value
        ShiftingCAimprovementThreshold = 0.05

        #Neither RebuildDefinition nor ElalXML will rebuild all the side chains
        SelectByOrigChainID = true
        Overlap = 0
eof
        if($status) then
          echo "   * Problem with SideAide. Cannot continue." | tee -a $LOG
          echo "COMMENT: error in chirality fixing" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                    >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          cd $BASE
          exit(1)
        else
          #Write debug messsages for building tools
          if ($LOCAL == 0 && `grep -c WRNG $WORKDIR/${PDBID}_chirfix.msg` > 0) then
            grep -H WRNG $WORKDIR/${PDBID}_chirfix.msg >> $DEBUGB
          endif
        endif

        #Reset the number of chirality fixes
        set NCHIRFX = `grep '#Residues to be rebuilt forcefully' $WORKDIR/$PDBID.todo | awk '{print $6}'`

        #If chiralities were fixed, consider model rebuilt
        if ($NCHIRFX > 0) then
          set BUILT = 1
        endif
      endif

      #Clean up
      rm -rf $WORKDIR/flip

    else
      #Give warning
      echo " o Validation-based rebuilding failed" | tee -a $LOG
      set WCERR = 1
      echo "COMMENT: WHAT_CHECK failed in rebuilding step." >> $DEBUG
      echo "PDB-REDO,$PDBID"                                >> $DEBUG
      cd $WORKDIR
      cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/${PDBID}_built.pdb
    endif
  else
    #Fill in some of the blanks
    set NSCFLIP = 0
    cp $WORKDIR/${PDBID}_scbuild.pdb $WORKDIR/${PDBID}_built.pdb
  endif

  #Present a rebuilding summary.
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Rebuilding details ******" | tee -a $LOG
  echo "Waters removed     : $NWATDEL" | tee -a $LOG
  echo "Peptides flipped   : $NBBFLIP" | tee -a $LOG
  echo "Side chains built  : $NSCBLT"  | tee -a $LOG
  echo "Side chains flipped: $NSCFLIP" | tee -a $LOG
  echo "Chirality fixes    : $NCHIRFX" | tee -a $LOG
endif

#Calculate the total number of chirality fixes
@ NCHIRFX = ($CHIFIX + $NCHIRFX)

############################################### Do another refinement run ################################################

echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Final refinement ******" | tee -a $LOG

#Skip the final refinement if the model was not rebuilt
unbuilt:

if ($BUILT == 0) then

  #Report
  echo "-Skipping the final refinement, because the model was not changed after re-refinement"  | tee -a $LOG

  #Copy some files
  cp $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/${PDBID}_final.pdb
  cp $WORKDIR/${PDBID}_besttls.mtz $WORKDIR/${PDBID}_final.mtz
  cp $WORKDIR/${PDBID}_besttls.log $WORKDIR/${PDBID}_final.log
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_refin.tls   $WORKDIR/${PDBID}_final.tls
  endif
  set BLTBEST = 'skip'

  #Copy refinement results
  set RFIN     = $RTLS
  set RFFIN    = $RFTLS
  set RFFINUNB = $RFTLSUNB
  set SIGRFFIN = $SIGRFTLS
  set RFFINZ   = $RFTLSZ

else
  #Update restraints
  #Nucleic acids (only if the previous run had no problems)
  if (-e $WORKDIR/nucleic.rest && $NPRUNEN == 0) then
    echo "-Updating nucleic acid restraints" | tee -a $LOG
    libg -p $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/stacking.rest -w sp >> $WORKDIR/libg.log
    libg -p $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/basepair.rest -w bp >> $WORKDIR/libg.log
    libg -p $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/chpucker.rest -w pu >> $WORKDIR/libg.log

    #Make Refmac use command files
    echo "external weight scale 5" >  $WORKDIR/nucleic.rest
    if (-e $WORKDIR/stacking.rest) then
      sed 's/  / /g' $WORKDIR/stacking.rest >> $WORKDIR/nucleic.rest
    endif
    if (-e $WORKDIR/basepair.rest) then
      sed 's/  / /g' $WORKDIR/basepair.rest >> $WORKDIR/nucleic.rest
    endif  
    if (-e $WORKDIR/chpucker.rest) then
      sed 's/  / /g' $WORKDIR/chpucker.rest >> $WORKDIR/nucleic.rest
    endif  
    echo "external weight scale 1" >> $WORKDIR/nucleic.rest
    set LIBGCMD = "@$WORKDIR/nucleic.rest"
  endif

  #Hydrogen bonds
  if (-e $WORKDIR/hbond.rest && ! -e $WORKDIR/homology.rest) then
    echo "-Updating hydrogen bond restraints" | tee -a $LOG

    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/${PDBID}_built.dssp >& $WORKDIR/dssp.log

    #Generate hydrogen bond restraints if a DSSP file can be made.
    if (! -e $WORKDIR/${PDBID}_built.dssp) then
      #DSSP failed. Cannot make restraints
      echo " o Cannot produce hydrogen bond restraints" | tee -a $LOG
      echo "DSSP: general error in second H-bond restraint generation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                                           >> $DEBUG
    else
      #Make the hydrogen bond restraints. PROGRAM: detectHbonds
      $TOOLS/detectHbonds -v \
      -pdb $WORKDIR/${PDBID}_built.pdb \
      -dssp $WORKDIR/${PDBID}_built.dssp \
      -output-name $PDBID \
      -tools $TOOLS >> $WORKDIR/hbondrest.log
      if ( ! -e $WORKDIR/${PDBID}_hbonds.rest) then
        echo " o Cannot generate hydrogen bond restraints for rebuilt model" | tee -a $LOG
        echo "detectHbonds: general error" >> $DEBUG
        echo "PDB-REDO,$PDBID"             >> $DEBUG
      else
        #Set up restraint commands
        sed 's/ [ ]*/ /g' $WORKDIR/${PDBID}_hbonds.rest > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
        set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
        set HBONDCMD    = "@$WORKDIR/hbond.rest"
        rm $WORKDIR/${PDBID}_built.dssp
      endif
    endif
  endif

  #Homology restraints
  if (-e $WORKDIR/homology.rest) then
    echo "-Updating homology-based restraints" | tee -a $LOG

    #Define the secondary structure
    mkdssp -i $WORKDIR/${PDBID}_built.pdb -o $WORKDIR/${PDBID}_built.dssp >& $WORKDIR/dssp.log

    #Make the restraints
    cd $WORKDIR/gbh

    $TOOLS/hoder -v \
    -alignout \
    -hitsummary \
    $HODERCMD \
    -pdb $WORKDIR/${PDBID}_built.pdb \
    -dssp $WORKDIR/${PDBID}_built.dssp \
    -blast $WORKDIR/$PDBID.blast \
    -fasta $WORKDIR/$PDBID.fasta \
    -output-name $PDBID \
    -output-dir $WORKDIR/gbh \
    -pdbdir $REDODIR \
    -edsdir $EDSDIR \
    -tools $TOOLS >> homologs.log

    #Set up the restraint commands
    #Homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/gbh/${PDBID}_homologyBased.restr > $WORKDIR/homology.rest #Replace one or more spaces by a single one
    set HOMOLWGTCMD = "EXTERNAL WEIGHT SCALE $HOMOLRESTRWGT"
    set HOMOLCMD    = "@$WORKDIR/homology.rest"

    #Non-homology hydrogen bonds
    sed 's/ [ ]*/ /g' $WORKDIR/gbh/${PDBID}_general.restr > $WORKDIR/hbond.rest #Replace one or more spaces by a single one
    set HBONDWGTCMD = "EXTERNAL WEIGHT SCALE $HBRESTRWGT"
    set HBONDCMD    = "@$WORKDIR/hbond.rest"

    #Go back down
    cd $WORKDIR
    rm $WORKDIR/${PDBID}_built.dssp
  endif

  #Metal restraints
  if ($DOMETALREST == 1 && `grep '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_built.pdb | cut -c 18-20 | grep -E -c 'ZN | ZN'` > 0 && $NPRUNEM == 0) then
    #Run Zen. PROGRAM: zen
    echo "-Updating metal restraints" | tee -a $LOG
    cp $WORKDIR/${PDBID}_built.pdb $WORKDIR/${PDBID}_built.old
    $TOOLS/zen \
    -v $WORKDIR/${PDBID}_built.old \
    $WORKDIR/${PDBID}_built.pdb \
    $WORKDIR/metal.raw \
    $WORKDIR/$PDBID.mtz \
    >> zen.log

    #Cleanup
    sed 's/  / /g' $WORKDIR/metal.raw > $WORKDIR/metal.rest
    rm $WORKDIR/metal.raw

    if (-e $WORKDIR/metal.rest) then
      set NMETALREST2 = `grep -c exte $WORKDIR/metal.rest`
      if ($NMETALREST2 != $NMETALREST) then
        echo " o Found different metal sites"      | tee -a $LOG
        echo " o Performing extra refinement using $NMETALREST2 restraints" | tee -a $LOG
      endif
    endif
  endif
  
  #Update antibody-specific restraints
  if (-e $WORKDIR/krab.rest) then
    echo "-Updating antibody-specific restraints" | tee -a $LOG
    #Source KRAB and create the virtual environment
    source $KRABSRC
      
    #Run KRAB
    krab --pyigclassify_root $PYIGCL $WORKDIR/${PDBID}_built.pdb
    
    #Filter with hbond and homology restraints if they exist
    if (-e $WORKDIR/homology.rest) then
      cp $WORKDIR/homology.rest $WORKDIR/homology.rest.unfiltered
      $TOOLS/filterest.py -v -o $WORKDIR/homology.rest $WORKDIR/krab.rest $WORKDIR/homology.rest.unfiltered >> $WORKDIR/filterest.log 
    endif
    if (-e $WORKDIR/hbond.rest) then
      cp $WORKDIR/hbond.rest $WORKDIR/hbond.rest.unfiltered
      $TOOLS/filterest.py -v -o $WORKDIR/hbond.rest $WORKDIR/krab.rest $WORKDIR/hbond.rest.unfiltered  >> $WORKDIR/filterest.log
    endif   
    
    #Switch off the KRAB virtual environment
    deactivate
  endif

  #set refinement parameters
  if ($TLSBEST == none) then
    #Autoweight
    set WGTTYPE = "AUTO"
    set WGRANGE = "2.50"
    set FCYCLE  = `echo $NCYCLE`
  else
    set FCYCLE  = 20
    #Increase the number of cycles if needed
    if ($NMETALREST2 != $NMETALREST) then
      @ FCYCLE = ($FCYCLE + 5)
    endif
    #Do short weight optimisation
    set WGTTYPE = "MATRIX"
    if ($TLSBEST == 5.00) then
      set WGRANGE = "5.00 7.00 3.00"
    else if ($TLSBEST == 3.00) then
      set WGRANGE = "3.00 4.00 2.00"
    else if ($TLSBEST == 2.00) then
      set WGRANGE = "2.00 2.50 1.50"
    else if ($TLSBEST == 1.50) then
      set WGRANGE = "1.50 1.80 1.00"
    else if ($TLSBEST == 1.00) then
      set WGRANGE = "1.00 1.30 0.70"
    else if ($TLSBEST == 0.70) then
      set WGRANGE = "0.70 0.80 0.50"
    else if ($TLSBEST == 0.50) then
      set WGRANGE = "0.50 0.60 0.30"
    else if ($TLSBEST == 0.30) then
      set WGRANGE = "0.30 0.40 0.10"
    else if ($TLSBEST == 0.10) then
      set WGRANGE = "0.10 0.15 0.05"
    else if ($TLSBEST == 0.05) then
      set WGRANGE = "0.05 0.07 0.03"
    else if ($TLSBEST == 0.03) then
      set WGRANGE = "0.03 0.05 0.02"
    else if ($TLSBEST == 0.01) then
      set WGRANGE = "0.01 0.02 .005"
    else if ($TLSBEST == .005) then
      set WGRANGE = ".005 .007 .002"
    else if ($TLSBEST == .002) then
      set WGRANGE = ".002 .005 .001"
    else if ($TLSBEST == .001) then
      set WGRANGE = ".001 .002 5e-4"
    else if ($TLSBEST == 5e-4) then
      set WGRANGE = "5e-4 .001 1e-4"
    else if ($TLSBEST == 1e-4) then
      set WGRANGE = "1e-4 2e-4 5e-5"
    else if ($TLSBEST == 1e-5) then
      set WGRANGE = "1e-5 2e-5 5e-6"
    else if ($TLSBEST == 1e-6) then
      set WGRANGE = "1e-6 2e-6 5e-7"
    else
      set WGRANGE = "1e-7 2e-7 5e-8"
    endif
  endif

  #Set up TLS and killswitch
  if ($DOTLS == 1 && $TLSUPDATE == 1) then
    set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_refin.tls TLSOUT $WORKDIR/${PDBID}_final.tls`
    #Do only a minor update of the TLS tensors (avoid instability)
    set TLSCMD     = `echo refi tlsc 5`
    set KILLSWITCH = "kill $TOOLS/pdb_redo.refmac"
  else
    set KILLSWITCH =    #Do not use the killswitch
    set TLSCMD     =    #Do not refine TLS
  endif

  #Run refmac with predefined matrix weights
  echo "-Final restraint weight optimisation with $WGTTYPE weights: $WGRANGE" | tee -a $LOG

  #Return label for running after TLS troubles
notlsblt:

  foreach WWGT (`echo $WGRANGE`)

    #Return label for job launching
bltrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_built.pdb \
      XYZOUT $WORKDIR/${PDBID}_blt$WWGT.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/${PDBID}_blt$WWGT.mtz \
      $LIBLIN \
      ${TLSFILS} \
      $SCATLIN \
<<eof >& $WORKDIR/${PDBID}_blt$WWGT.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        $TLSCMD
        tlsd waters exclude
        ncyc $FCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $WGTSIG $WGTTYPE $WWGT
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $OCCCMD
        $HARMCMD
        $METALCMD
        $LIBGCMD
        $RESTCMD
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES $ANOMCOEF
        $ANOMCMD
        pdbout copy remarks 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KRABWGTCMD
        $KRABCMD
        $KILLSWITCH
        END
eof
    else
      #Wait a bit to start again
      sleep 10
      goto bltrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Check for problems by seeing if the output mtz file exists
  foreach WWGT (`echo $WGRANGE`)
    if (! -e $WORKDIR/${PDBID}_blt$WWGT.mtz) then
      #Was the killswitch activated
      if (`grep -a -c 'Program terminated by user' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
        #Yes, check if TLS was used.
        if ($DOTLS == 1 && $TLSUPDATE == 1) then
          #Try running without updating the TLS tensors.
          set TLSFILS    = `echo TLSIN $WORKDIR/${PDBID}_refin.tls`
          set TLSCMD     =
          set TLSUPDATE  = 0
          set KILLSWITCH =    #Do not use the killswitch

          #Give warning and write debug statement.
          echo " " | tee -a $LOG
          echo " o The refinement was unstable." | tee -a $LOG
          echo "   * Trying refinement without updating the TLS model " | tee -a $LOG

          #Rerun Refmac without updating the TLS
          goto notlsblt
        else if
          #Check for atom naming problems
          if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
            #Is this caused by an alternative compound
            set RESN = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 39-42`
            set CHID = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 45-45`
            #Check whether there is more than 1 residue name for a given chain and residue number.
            if (`grep -E $CHID' {0,3}'$RESN $WORKDIR/${PDBID}_besttls.pdb | cut -c 18-20 | sort -u | wc -l` > 1) then
             
              #Write warning
              echo " " | tee -a $LOG
              echo " o The rebuilt model caused problems in Refmac" | tee -a $LOG
              echo "   * Falling back to the re-refined model"      | tee -a $LOG
              echo "COMMENT: refmac: could not use rebuilt model" >> $DEBUG
              echo "PDB-REDO,$PDBID"                              >> $DEBUG
              
              #Undo the rebuilding
              set NWATDEL = 0
              set NBBFLIP = 0
              set NSCBLT  = 0
              set NSCFLIP = 0
              set NCHIRFX = $CHIFIX
              set BUILT   = 0

              #Start again form the pre-rebuilding model
              goto unbuilt
              
            else
              #Fatal atom naming problem
              echo " " | tee -a $LOG
              echo " o Problem with restraint generation in refmac. Cannot continue." | tee -a $LOG
              echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
              echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
              if ($SERVER == 1) then
                #Write out status files
                touch $STDIR/stoppingProcess.txt
                touch $STDIR/processStopped.txt
              endif
              exit(1)
            endif
          endif
        else
          #TLS was not the problem. The error is fatal.
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          exit(1)
        endif
      else
        #Check for atom naming problems
        if (`grep -a -c 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log` != 0) then
          #Is this caused by an alternative compound
          set RESN = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 39-42`
          set CHID = `grep -a 'is absent in the library' $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | cut -c 45-45`
          #Check whether there is more than 1 residue name for a given chain and residue number.
          if (`grep -E $CHID' {0,3}'$RESN $WORKDIR/${PDBID}_besttls.pdb | cut -c 18-20 | sort -u | wc -l` > 1) then
           
            #Write warning
            echo " " | tee -a $LOG
            echo " o The rebuilt model caused problems in Refmac" | tee -a $LOG
            echo "   * Falling back to the re-refined model"      | tee -a $LOG
            echo "COMMENT: refmac: could not use rebuilt model" >> $DEBUG
            echo "PDB-REDO,$PDBID"                              >> $DEBUG
              
            #Undo the rebuilding
            set NWATDEL = 0
            set NBBFLIP = 0
            set NSCBLT  = 0
            set NSCFLIP = 0
            set NCHIRFX = $CHIFIX
            set BUILT   = 0

            #Start again form the pre-rebuilding model
            goto unbuilt
        
          else
            #Fatal atom naming problem
            echo " " | tee -a $LOG
            echo " o Problem with restraint generation in refmac. Cannot continue." | tee -a $LOG
            echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
            echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
            if ($SERVER == 1) then
              #Write out status files
              touch $STDIR/stoppingProcess.txt
              touch $STDIR/processStopped.txt
            endif
            exit(1)
          endif
        else
          #Something else is wrong, give an error mesage and quit.
          echo " " | tee -a $LOG
          echo " o Problem with refmac. Cannot continue." | tee -a $LOG
          echo "COMMENT: refmac: error in post-rebuilding refinement" >> $WHYNOT
          echo "PDB-REDO,$PDBID"                                      >> $WHYNOT
          if ($SERVER == 1) then
            #Write out status files
            touch $STDIR/stoppingProcess.txt
            touch $STDIR/processStopped.txt
          endif
          cd $BASE
          exit(1)
        endif  
      endif
    endif
  end

  #Report
  foreach WWGT (`echo $WGRANGE`)
    set TFREE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1 | awk '{print $3}'`
    echo " o Performed ${FCYCLE}-cycle refinement with $WGTTYPE weight $WWGT (R-free = $TFREE)" | tee -a $LOG
  end

  #Give refinement output
  #No known errors in the refinement
  set BLTERR  = 0

  #If the refinement was done with automated weighting, just rename the files...
  if ($WGTTYPE == "AUTO") then
    set BLTBEST = 'auto'
    mv $WORKDIR/${PDBID}_blt$WWGT.pdb $WORKDIR/${PDBID}_final.pdb
    mv $WORKDIR/${PDBID}_blt$WWGT.mtz $WORKDIR/${PDBID}_final.mtz
    mv $WORKDIR/${PDBID}_blt$WWGT.log $WORKDIR/${PDBID}_final.log
  else
    #Create second.refi file
    if ($GOTR == 0 || $ZCALERR == 1) then
      echo "$RCAL $RFCALUNB" > $WORKDIR/${PDBID}blt.refi  #Use recalculated R and 'unbiased' R-free as benchmarks
    else
      echo "$RCAL $RFCAL" > $WORKDIR/${PDBID}blt.refi  #Use R and R-free obtained from the recalculation as benchmarks
    endif
    foreach WWGT (`echo $WGRANGE`)
      if (-e $WORKDIR/${PDBID}_blt$WWGT.pdb) then
        set LINE = `tail -n $LOGSTEP $WORKDIR/${PDBID}_blt$WWGT.log | head -n 1`
        echo "$WWGT $LINE" >> $WORKDIR/${PDBID}blt.refi
      else
        echo " o $WORKDIR/${PDBID}_blt$WWGT.pdb is missing" | tee -a $LOG
        set BLTERR = 1
      endif
    end

    #Pick the best re-refined structure
    set BLTBEST = `$TOOLS/picker -s -f $ESTRICT $WORKDIR/${PDBID}blt.refi $NTSTCNT $RFRRAT $RMSZB $RMSZA`
    cp $WORKDIR/${PDBID}_blt$BLTBEST.pdb $WORKDIR/${PDBID}_final.pdb
    cp $WORKDIR/${PDBID}_blt$BLTBEST.mtz $WORKDIR/${PDBID}_final.mtz
    cp $WORKDIR/${PDBID}_blt$BLTBEST.log $WORKDIR/${PDBID}_final.log
  endif

  #Set R(-free)
  set RFIN  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_final.log | head -n 1 | awk '{print $2}'`
  set RFFIN = `tail -n $LOGSTEP $WORKDIR/${PDBID}_final.log | head -n 1 | awk '{print $3}'`
  setenv PRFIN  $RFIN   #Set value for bits of Perl code in the script
  setenv PRFFIN $RFFIN  #Set value for bits of Perl code in the script

  #Validate R-free
  set SIGRFFIN = `perl -e 'printf ("%.4f\n", $ENV{PRFFIN}/sqrt($ENV{PNTSTCNT}));'`
  setenv PSIGRFFIN $SIGRFFIN
  set RFFINUNB = `perl -e 'printf ("%.4f\n", $ENV{PRFIN}*$ENV{PRFRRAT});'`
  setenv PRFFINUNB $RFFINUNB
  set RFFINZ   = `perl -e 'printf ("%.2f\n", ($ENV{PRFFINUNB}-$ENV{PRFFIN})/$ENV{PSIGRFFIN});'`

  echo " " | tee -a $LOG
  echo " " | tee -a $LOG

  echo "****** Final refinement details ******" | tee -a $LOG
  if ($WGTTYPE == 'MATRIX') then
    echo "Best matrix weight: $BLTBEST" | tee -a $LOG
  endif
  echo "Resulting R-factor: $RFIN"     | tee -a $LOG
  echo "Resulting R-free  : $RFFIN"    | tee -a $LOG
  echo "'Unbiased' R-free : $RFFINUNB" | tee -a $LOG
  echo "sigma(R-free)     : $SIGRFFIN" | tee -a $LOG
  echo "R-free Z-score    : $RFFINZ"   | tee -a $LOG
endif

#Copy back the SEQRES records
if ($GOTSEQRES == 1) then
  cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final.bak
  #Run seqrescopier. PROGRAM: seqrescopier
  $TOOLS/seqrescopier -v \
  -pdbinw  $WORKDIR/pdb$PDBID.ent \
  -pdbinwo $WORKDIR/${PDBID}_final.bak \
  -pdbout  $WORKDIR/${PDBID}_final.pdb > seqrescopier.log 
else
  #Do nothing
endif 

#Run flipper again
cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final.bak
$TOOLS/flipper -v \
-pdbin  $WORKDIR/${PDBID}_final.bak \
-pdbout $WORKDIR/${PDBID}_final.pdb >> $WORKDIR/flipper.log

###########################################  Correlation coefficients  ###################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG

#Extract the data
set CCWOLD = `grep 'CORRELATION COEFFICIENT FO-FC     ' $WORKDIR/${PDBID}_0cyc.pdb  | awk '{print $7}'`
set CCWFIN = `grep 'CORRELATION COEFFICIENT FO-FC     ' $WORKDIR/${PDBID}_final.pdb | awk '{print $7}'`
set CCFOLD = `grep 'CORRELATION COEFFICIENT FO-FC FREE' $WORKDIR/${PDBID}_0cyc.pdb  | awk '{print $8}'`
set CCFFIN = `grep 'CORRELATION COEFFICIENT FO-FC FREE' $WORKDIR/${PDBID}_final.pdb | awk '{print $8}'`

#Get the Z-scores
set ZCCW = `echo "$CCWOLD $CCWFIN $PWORKCNT" | awk '{zcco = 0.5*(log((1+$1)/(1-$1))); zccf = 0.5*(log((1+$2)/(1-$2))); zchange = (zccf-zcco)/sqrt(2/($3-3))} END {printf "%6.2f\n", zchange}'`
set ZCCF = `echo "$CCFOLD $CCFFIN $PNTSTCNT" | awk '{zcco = 0.5*(log((1+$1)/(1-$1))); zccf = 0.5*(log((1+$2)/(1-$2))); zchange = (zccf-zcco)/sqrt(2/($3-3))} END {printf "%6.2f\n", zchange}'`

#Report the results
echo "****** Reciprocal space correlation ******" | tee -a $LOG
echo '                              Before  Final  Z-score'    | tee -a $LOG
echo "Work correlation coefficient: $CCWOLD   $CCWFIN  $ZCCW"  | tee -a $LOG
echo "Free correlation coefficient: $CCFOLD   $CCFFIN  $ZCCF"  | tee -a $LOG



############################################   Full cross validation  ####################################################

#Do full cross-validation if the R-free set is too small or if asked by the user...
# ... but only if there was a successful refinement
if ( ($NTSTCNT < 500 || $CROSS == 1) && !($BLTBEST == 'skip' && $TLSBEST == 'none')) then
  echo " " | tee -a $LOG
  echo " " | tee -a $LOG
  echo "****** Cross validation ******" | tee -a $LOG
  if ($NTSTCNT < 500) then
    echo "-The test set contained fewer than 500 reflections" | tee -a $LOG
    set CROSS = 1
  endif

  #Calculate K
  set KFOLD = `mtzdmp $WORKDIR/$PDBID.mtz -n 0 | grep '  FREE' | cut -c 21-22 | awk '{print $1+1}'`
  @ MAXSET  = ($KFOLD - 1)
  echo "-Performing $KFOLD-fold cross validation" | tee -a $LOG

  #Set up the refinement
  if ($DOTLS == 1) then
    set TLSFILS = `echo TLSIN $WORKDIR/${PDBID}_final.tls`
  endif
  set BCMD = "bfac set $BSET"
  @ NCYCLE = ($NCYCLE + 40)
  #Keep the final refinement settings unless there was no final refinement
  if ($BLTBEST == 'skip') then
    set CWEIGHT = `echo "$WGTSIG MATRIX $TLSBEST"`
  else if ($BLTBEST == 'auto') then
    set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
  else
    set CWEIGHT = `echo "$WGTSIG $WGTTYPE $BLTBEST"`
  endif

  #Perturb the atomic coordinates if the B-factor model is OVER
  if ($BREFTYPE == "OVER") then
    #Run pdbset
    pdbset \
    XYZIN $WORKDIR/${PDBID}_final.pdb \
    XYZOUT $WORKDIR/${PDBID}_crossin.pdb \
    <<eof > $WORKDIR/pdbset.log
      NOISE 0.05
eof
  else
    cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_crossin.pdb
  endif

  #Do the refinements
  foreach SET (`seq -s " " 0 $MAXSET`)

    #Return label for job launching
xvalrunning:

    #Only launch new jobs when the number of cores is not exceeded
    #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
    jobs > $WORKDIR/jobs.log

    if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

      echo " o Starting refinement with test set $SET" | tee -a $LOG

      refmac5 \
      XYZIN  $WORKDIR/${PDBID}_crossin.pdb \
      XYZOUT $WORKDIR/cross$SET.pdb \
      HKLIN  $WORKDIR/$PDBID.mtz \
      HKLOUT $WORKDIR/cross$SET.mtz \
      $LIBLIN \
      ${TLSFILS} \
      $SCATLIN \
  <<eof > $WORKDIR/${PDBID}_cross$SET.log &
        $SCATTERCMD
        make check NONE
        make hydrogen $HYDROGEN hout NO peptide NO cispeptide YES -
          ssbridge $SSBOND $RSYMM $SUGAR $CONNECTIVITY link NO
        refi type REST resi MLKF meth CGMAT bref $BREFTYPE
        $REFIRES
        free $SET
        $BCMD
        tlsd waters exclude
        ncyc $NCYCLE
        scal type $SOLVENT $SCALING
        solvent YES
        $MASKPAR
        $LOWMEM
        weight $CWEIGHT
        monitor MEDIUM torsion 10.0 distance 10.0 angle 10.0 plane 10.0 -
          chiral 10.0  bfactor 10.0 bsphere  10.0 rbond 10.0 ncsr  10.0
        $TWIN
        $NCSTYPE
        $NCSALIGN
        $NCSNEIGH
        $NCSSTRICT
        $JELLY
        $TORSION
        $HARMCMD
        $METALCMD
        $LIBGCMD
        $RESTCMD
        temp $BBEST
        blim 2.0 999.0
        labin FP=FP SIGFP=SIGFP FREE=FREE $PHASES
        pdbout copy remarks 280 350
        NOHARVEST
        $HOMOLWGTCMD
        $HOMOLCMD
        $HBONDWGTCMD
        $HBONDCMD
        $KRABWGTCMD
        $KRABCMD
        END
eof
     else
      #Wait a bit to start again
      sleep 10
      goto xvalrunning
    endif
  end

  #Wait for the jobs to finish
  wait

  #Make a summary only if the refinement was successful
  echo "-Analysing cross validation results" | tee -a $LOG
  foreach SET (`seq -s " " 0 $MAXSET`)
    if (! -e $WORKDIR/cross$SET.mtz) then
      echo " o Problem in cross validation using test set $SET" | tee -a $LOG
      echo "COMMENT: refmac: error in cross validation" >> $DEBUG
      echo "PDB-REDO,$PDBID"                            >> $DEBUG
    else
      set LINE  = `tail -n $LOGSTEP $WORKDIR/${PDBID}_cross$SET.log | head -n 1`
      set NTREF = `mtzdmp $WORKDIR/$PDBID.mtz -n -1 | grep -v '?' | grep -E "^.{4}[0123456789].{3}[0123456789]" | awk '{printf "%d\n", $4}' | grep -c "$SET"`
      echo "$SET $LINE $NTREF" >> $WORKDIR/${PDBID}.rtest
    endif
  end

  #Analyse results and show them

  #R-complete
  foreach SET (`seq -s " " 0 $MAXSET`)

    #Copy back the testset flags (circumvents a Refmac quirk)
    cad \
    HKLIN1 $WORKDIR/cross${SET}.mtz \
    HKLIN2 $WORKDIR/$PDBID.mtz \
    HKLOUT $WORKDIR/crossf${SET}.mtz \
    <<eof >> $WORKDIR/rcomplete.log
    LABIN FILE 1  E1=FP E2=FC_ALL_LS E3=PHIC_ALL_LS
    LABIN FILE 2  E1=FREE
    END
eof

    #Convert the file to mmCIF
    mtz2various \
    HKLIN  $WORKDIR/crossf${SET}.mtz \
    HKLOUT $WORKDIR/cross${SET}.cif << eof >> $WORKDIR/rcomplete.log
    OUTPUT CIF data_rcom
    FREEVAL $SET.0
    LABIN FP=FP FC=FC_ALL_LS PHIC=PHIC_ALL_LS FREE=FREE
    END
eof

    #Extract and consolidate the data
    grep ' 1 1 1' $WORKDIR/cross${SET}.cif | grep ' f ' >> $WORKDIR/rcomplete.hkl
  end

  #Get the R-complete
  set RCOMPLETE =  `cat $WORKDIR/rcomplete.hkl | awk '{if ($8 > $9) {sumdif += ($8 - $9)} else {sumdif += ($9 - $8)}; sumf += $8} END {printf "%6.4f", sumdif/sumf}'`
  echo " o R-complete = $RCOMPLETE" | tee -a $LOG

  #Averages
  $TOOLS/longinus -v \
  $WORKDIR/$PDBID.rtest > $WORKDIR/longinus.log

  #Did the crossvalidation work
  if (`grep -c 'Could not find any valid refinement data' $WORKDIR/longinus.log` != 0) then
    #Cross validation failed write an error message
    set CROSS = 0
    echo "-Cross validation failed" | tee -a $LOG
    echo "COMMENT: error in cross validation " >> $DEBUG
    echo "PDB-REDO,$PDBID"                     >> $DEBUG
  else
    #Crossvalidation successful. Write out details.
    echo " " | tee -a $LOG
    echo "****** Cross validation details ******" | tee -a $LOG
    tail -n 3 $WORKDIR/longinus.log | tee -a $LOG

    #Mine the values
    set CRFACT  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 10-15`
    set CSRFACT = `tail -n 1 $WORKDIR/longinus.log | cut -c 10-15`
    set CRFREE  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 20-25`
    set CSRFREE = `tail -n 1 $WORKDIR/longinus.log | cut -c 20-25`
    set CGAP    = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 28-33`
    set CSGAP   = `tail -n 1 $WORKDIR/longinus.log | cut -c 28-33`
    set CNTEST  = `tail -n 2 $WORKDIR/longinus.log | head -n 1 | cut -c 35-40`
    set CSNTEST = `tail -n 1 $WORKDIR/longinus.log | cut -c 35-40`

    #Clean up
    foreach SET (`seq -s " " 0 $MAXSET`)
      rm $WORKDIR/cross$SET.mtz
      rm $WORKDIR/crossf$SET.mtz
      rm $WORKDIR/cross$SET.pdb
      rm $WORKDIR/cross$SET.cif
    end
  endif
endif

############################################   Validate structure   ######################################################
echo " " | tee -a $LOG
echo " " | tee -a $LOG
echo "****** Validation details ******" | tee -a $LOG
echo "-Validating structure" | tee -a $LOG

#Initialise
set WCERR = 0

#Go to temporary running directory
setenv WCWORF $WORKDIR/wctemf
mkdir -p $WCWORF
cd $WCWORF

#Get the PDB file
cp $WORKDIR/${PDBID}_final.pdb $WCWORF

#Do the actual validation
$WC/bin/whatcheck $WCWORF/${PDBID}_final.pdb Y Y Y >& $WORKDIR/wc_final.log

#Check for an output file
if (-e $WCWORF/pdbout.txt) then
  #Do Nothing
else
  #Give warning
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate final structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                     >> $DEBUG
endif

#Create webpage
$WC/bin/pdbout2html >>& $WORKDIR/wc_final.log

#Check index.html completeness
if (`grep -c "Final summary" $WCWORF/pdbout.html` == 0 && `grep -c "Summary report" $WCWORF/pdbout.html` == 0) then
  echo " o Validation failed" | tee -a $LOG
  set WCERR = 1
  echo "COMMENT: WHAT_CHECK cannot validate final structure" >> $DEBUG
  echo "PDB-REDO,$PDBID"                                     >> $DEBUG
endif

#Clean up on aisle six
mkdir -p $WORKDIR/wf

if (-e $WCWORF/pdbout.txt) then
  mv $WCWORF/pdbout.txt $WORKDIR/wf/
endif
if (-e $WCWORF/check.db) then
  bzip2 -f $WCWORF/check.db
  mv $WCWORF/check.db.bz2 $WORKDIR/wf/
endif
#server only
if ($SERVER == 1) then
  mv $WCWORF/pdbout.html $WORKDIR/wf/index.html
  mv $WCWORF/*.gif $WORKDIR/wf/ >& /dev/null
endif

#Go to start directory
cd $WORKDIR
rm -rf $WCWORF

#Extract stats from refitted structure
set FNATOM = `grep -c -E '^[AH][TE][OT][MA]' $WORKDIR/${PDBID}_final.pdb`
if (-e $WORKDIR/wf/pdbout.txt) then
  if (`grep -a -c '1st generation packing quality :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZPAK1 = 'NA'
  else
    set FZPAK1 = `grep -a '1st generation packing quality :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c '2nd generation packing quality :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZPAK2 = 'NA'
  else
    set FZPAK2 = `grep -a '2nd generation packing quality :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Ramachandran plot appearance   :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FZRAMA = 'NA'
  else
    set FZRAMA = `grep -a 'Ramachandran plot appearance   :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FCHI12 = 'NA'
  else
    set FCHI12 = `grep -a 'chi-1/chi-2 rotamer normality  :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Backbone conformation          :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FBCONF = 'NA'
  else
    set FBCONF = `grep -a 'Backbone conformation          :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond lengths                   :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FBRMSZ = 'NA'
  else
    set FBRMSZ = `grep -a 'Bond lengths                   :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  if (`grep -a -c 'Bond angles                    :' $WORKDIR/wf/pdbout.txt` == 0) then
    set FARMSZ = 'NA'
  else
    set FARMSZ = `grep -a 'Bond angles                    :' $WORKDIR/wf/pdbout.txt | head -n 1 | cut -c 36-42`
  endif
  #Get bump scores
  if ($OLDWC == 1) then
    set FWBMPS = 'NA'
    if (`grep -a -A1 -E '^.{31,32}<-->' $WORKDIR/wf/pdbout.txt | tail -n 1 | grep -c 'total of'` == 0) then
      set FBUMPS = `grep -a -c -E '^.{31,32}<-->' $WORKDIR/wf/pdbout.txt`
    else
      set FBUMPS = `grep -a -A1 -E '^.{31,32}<-->' $WORKDIR/wf/pdbout.txt | tail -n 1 | awk '{print $8}'`
    endif
  else
    if (`grep -a -c 'Total number of bumps:' $WORKDIR/wf/pdbout.txt` == 0) then
      set FBUMPS = 0
      set FSBMPL = 0
      set FWBMPS = 0.000
    else
      set FBUMPS = `grep -a 'Total number of bumps:' $WORKDIR/wf/pdbout.txt | cut -c 24-28`
      set FSBMPL = `grep -a 'Total squared bump value:' $WORKDIR/wf/pdbout.txt | cut -c 27-33`
      set FWBMPS = `echo $FNATOM $FSBMPL | awk '{printf "%5.3f\n", 100*$2/$1}'`
    endif
  endif
  #Count unsatisfied donors and acceptors
  if (`grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set FHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set FHBDON = `grep -a -A 80 'The buried hydrogen bond donors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  if (`grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c 'total of'` == 0) then
    set FHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep -c -E '^.{11}\('`
  else
    set FHBACC = `grep -a -A 80 'The buried side-chain hydrogen bond acceptors' $WORKDIR/wf/pdbout.txt | grep -B 80 -m 1 '#' | grep 'total of' | awk '{print $8}'`
  endif
  @ FHBUNS = ($FHBDON + $FHBACC)
  if (`grep -a -c 'Buried donors:' $WORKDIR/wf/pdbout.txt` == 0) then
    set FHBSAT = 'NA'
  else
    @ FBHBDA = (`grep -a 'Buried donors:' $WORKDIR/wf/pdbout.txt | awk '{print $3}'` + `grep -a 'Buried acceptors:' $WORKDIR/wf/pdbout.txt | awk '{print $3}'`)
    set FHBSA1 = `grep -a 'with a H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $5} END {print SUM}'`
    set FHBSA2 = `grep -a 'with a poor H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $6} END {print 0.5*SUM}'`
    set FHBSA3 = `grep -a 'with only a very poor H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $8} END {print 0.25*SUM}'`
    set FHBSA4 = `grep -a 'essentially without H-bond:' $WORKDIR/wf/pdbout.txt | awk '{SUM += $5} END {print 0.125*SUM}'`
    set FHBSAT = `echo $FHBSA1 $FHBSA2 $FHBSA3 $FHBSA4 $FBHBDA | awk '{printf "%5.3f\n" ,($1 + $2 + $3 + $4)/$5}'`
  endif

else
  set FZPAK1 = 'NA'
  set FZPAK2 = 'NA'
  set FZRAMA = 'NA'
  set FCHI12 = 'NA'
  set FBCONF = 'NA'
  set FBRMSZ = 'NA'
  set FARMSZ = 'NA'
  set FBUMPS = 'NA'
  set FWBMPS = 'NA'
  set FHBUNS = 'NA'
  set FHBSAT = 'NA'
endif

#Calculate percentiles vs the PDB
if ($FZRAMA == 'NA') then
  set TFZRAMA = 'NA'
else
  set TFZRAMA = `cat $TOOLS/zrama.sort | awk -v VAL=$FZRAMA '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZRAMA == 'NA') then
  set TOZRAMA = 'NA'
else
  set TOZRAMA = `cat $TOOLS/zrama.sort | awk -v VAL=$OZRAMA '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FCHI12 == 'NA') then
  set TFCHI12 = 'NA'
else
  set TFCHI12 = `cat $TOOLS/chi12.sort | awk -v VAL=$FCHI12 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OCHI12 == 'NA') then
  set TOCHI12 = 'NA'
else
  set TOCHI12 = `cat $TOOLS/chi12.sort | awk -v VAL=$OCHI12 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FWBMPS == 'NA') then
  set TFWBMPS = 'NA'
else
  set TFWBMPS = `cat $TOOLS/bumps.sort | awk -v VAL=$FWBMPS '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OWBMPS == 'NA') then
  set TOWBMPS = 'NA'
else
  set TOWBMPS = `cat $TOOLS/bumps.sort | awk -v VAL=$OWBMPS '{if ($1 > VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FZPAK2 == 'NA') then
  set TFZPAK2 = 'NA'
else
  set TFZPAK2 = `cat $TOOLS/zpak2.sort | awk -v VAL=$FZPAK2 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZPAK2 == 'NA') then
  set TOZPAK2 = 'NA'
else
  set TOZPAK2 = `cat $TOOLS/zpak2.sort | awk -v VAL=$OZPAK2 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FZPAK1 == 'NA') then
  set TFZPAK1 = 'NA'
else
  set TFZPAK1 = `cat $TOOLS/zpak1.sort | awk -v VAL=$FZPAK1 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OZPAK1 == 'NA') then
  set TOZPAK1 = 'NA'
else
  set TOZPAK1 = `cat $TOOLS/zpak1.sort | awk -v VAL=$OZPAK1 '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($FHBSAT == 'NA') then
  set TFHBSAT = 'NA'
else
  set TFHBSAT = `cat $TOOLS/hbsat.sort | awk -v VAL=$FHBSAT '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif
if ($OHBSAT == 'NA') then
  set TOHBSAT = 'NA'
else
  set TOHBSAT = `cat $TOOLS/hbsat.sort | awk -v VAL=$OHBSAT '{if ($1 < VAL) RANK = NR} END {printf "%3.0f", RANK*100/NR}'`
endif

#Run FoldX (if available)
if ($?FOLDX) then
  echo "-Running FoldX" | tee -a $LOG
  #Run the final PDB file
  $FOLDX --command=Stability --pdb-dir=$WORKDIR --pdb=${PDBID}_final.pdb --output-file=$PDBID >>& $WORKDIR/foldx.log

  #Set the value
  if (`grep -c 'Total          =' $WORKDIR/foldx.log` != 0) then
    set FGFOLD = `grep 'Total          =' $WORKDIR/foldx.log | tail -n 1 | awk '{print $3}'`
  else
    set FGFOLD = 'NA'
  endif
else
  set FGFOLD = 'NA'
endif

#Run distel if needed (note that the restraints may have changed so the rmsZ values may be different)
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  #Consolidate the distance restraints
  if (-e $WORKDIR/homology.rest) then
    cat $WORKDIR/homology.rest > $WORKDIR/allhb.rest
  endif
  if (-e $WORKDIR/hbond.rest) then
    cat $WORKDIR/hbond.rest >> $WORKDIR/allhb.rest
  endif
  
  #Calculate rmsZ values. Also make YASARA scene file for final model.
  set OHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set NHRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  set FHRMSZ = `$TOOLS/distel.py -m $WORKDIR/${PDBID}_final.pdb $WORKDIR/allhb.rest | tail -n 1 | awk '{print $4}'`
  cp $WORKDIR/distel_distance_restraints.mcr $WORKDIR/hbond.mcr
else 
  set OHRMSZ = 'NA'
  set NHRMSZ = 'NA'
  set FHRMSZ = 'NA'
endif

#Run distel for torsion angle restraints
if ($NKRABREST > 0) then
  set OKRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_0cyc.pdb $WORKDIR/krab.rest | tail -n 1 | awk '{print $5}'`
  set NKRMSZ = `$TOOLS/distel.py $WORKDIR/${PDBID}_besttls.pdb $WORKDIR/krab.rest | tail -n 1 | awk '{print $5}'`
  set FKRMSZ = `$TOOLS/distel.py -m $WORKDIR/${PDBID}_final.pdb $WORKDIR/krab.rest | tail -n 1 | awk '{print $5}'`
  cp $WORKDIR/distel_torsion_restraints.mcr $WORKDIR/krab.mcr
else
  set OKRMSZ = 'NA'
  set NKRMSZ = 'NA'
  set FKRMSZ = 'NA'
endif

#Print results
echo '                                     Before  Re-ref  Final' | tee -a $LOG
echo "1st generation packing quality     : $OZPAK1 $NZPAK1 $FZPAK1" | tee -a $LOG
echo "2nd generation packing quality     : $OZPAK2 $NZPAK2 $FZPAK2" | tee -a $LOG
echo "Ramachandran plot appearance       : $OZRAMA $NZRAMA $FZRAMA" | tee -a $LOG
echo "chi-1/chi-2 rotamer normality      : $OCHI12 $NCHI12 $FCHI12" | tee -a $LOG
echo "Backbone conformation              : $OBCONF $NBCONF $FBCONF" | tee -a $LOG
echo " " | tee -a $LOG
echo "Bond length RMS Z-score            : $OBRMSZ $NBRMSZ $FBRMSZ" | tee -a $LOG
echo "Bond angle RMS Z-score             : $OARMSZ $NARMSZ $FARMSZ" | tee -a $LOG
echo " " | tee -a $LOG
echo "Total number of bumps              : $OBUMPS $NBUMPS $FBUMPS" | tee -a $LOG
echo "Weighted bump severity score       : $OWBMPS $NWBMPS $FWBMPS" | tee -a $LOG
echo " " | tee -a $LOG
echo "Unsatisfied H-bond donors/acceptors: $OHBUNS $NHBUNS $FHBUNS" | tee -a $LOG
echo "H-bond satisfaction fraction       : $OHBSAT $NHBSAT $FHBSAT" | tee -a $LOG
if ($NHOMOLREST > 0 || $NHBONDREST > 0) then
  echo "H-bond restraint RMS Z-score       : $OHRMSZ $NHRMSZ $FHRMSZ" | tee -a $LOG
endif
if ($NKRABREST > 0) then
  echo "KRAB restraint RMS Z-score         : $OKRMSZ $NKRMSZ $FKRMSZ" | tee -a $LOG
endif
if ($?FOLDX) then
  echo " " | tee -a $LOG
  echo "Gibbs folding energy (kcal/mol)    : $OGFOLD $NGFOLD $FGFOLD" | tee -a $LOG
endif
echo " " | tee -a $LOG

#Give the side chain details
echo " " | tee -a $LOG
echo "****** Protein side chain details ******" | tee -a $LOG

#Set fallback values
set NDROTA = 'NA'
set HBFLIP = 'NA'
set STFLIP = 'NA'

#Get rotamer and peptide state changes
if ($GOT_PROT == T) then
  #Check for changed rotamers
  echo "-Analysing torsion angles with rotacompare" | tee -a $LOG
  $TOOLS/rotacompare -v \
  -pdb1 $WORKDIR/pdb${PDBID}.ent \
  -pdb2 $WORKDIR/${PDBID}_final.pdb \
  -output-name $PDBID \
  -output-dir $WORKDIR >> $WORKDIR/${PDBID}_torsions.log
  if (-e $WORKDIR/${PDBID}_torsiondiff.txt) then
    #Rotacompare worked
  else
    #Give error message
    echo " o Problem with rotacompare" | tee -a $LOG
    echo "COMMENT: rotacompare: error in rotamer validation" >> $DEBUG
    echo "PDB-REDO,$PDBID"                                   >> $DEBUG

    #Make and empty file
    touch $WORKDIR/${PDBID}_torsiondiff.txt
  endif
endif

#Create the coot_tour
$TOOLS/coot_tour -v $WORKDIR/pdb$PDBID.ent $WORKDIR/${PDBID}_final.pdb \
  $WORKDIR/${PDBID}_final.scm $WORKDIR/${PDBID}_final.py \
  $WORKDIR/${PDBID}_pepflip.log $WORKDIR/${PDBID}_torsiondiff.txt  > $WORKDIR/coot_tour.log

#Count things
set HBFLIP = `grep 'H-bond flips marked' $WORKDIR/coot_tour.log | awk '{print $1}'`
set NDROTA = `grep 'changed rotamers marked' $WORKDIR/coot_tour.log | awk '{print $1}'`
set CTFLIP = `grep 'cis-trans isomerisations marked' $WORKDIR/coot_tour.log | awk '{print $1}'`
set TCFLIP = `grep 'trans-cis isomerisations marked' $WORKDIR/coot_tour.log | awk '{print $1}'`
set PEPFIX = `grep 'removed distorted peptides marked' $WORKDIR/coot_tour.log | awk '{print $1}'`
set PEPDIS = `grep 'introduced distorted peptides marked' $WORKDIR/coot_tour.log | awk '{print $1}'`

#Give the summary
echo " " | tee -a $LOG
echo "Changed rotamers     : $NDROTA" | tee -a $LOG
echo "Hydrogen bond flips  : $HBFLIP" | tee -a $LOG
echo "Cis-trans flips      : $CTFLIP" | tee -a $LOG
echo "Trans-cis flips      : $TCFLIP" | tee -a $LOG
echo "Peptides fixed       : $PEPFIX" | tee -a $LOG
echo "Peptides distorted   : $PEPDIS" | tee -a $LOG
echo " " | tee -a $LOG



#Give the density fit details
echo " " | tee -a $LOG
echo "****** Electron density map details ******" | tee -a $LOG

#Set fallback values
set RSRB  = 'NA'
set RSRW  = 'NA'
set RSCCB = 'NA'
set RSCCW = 'NA'
set PLOTS = 0


echo "-Analysing density maps with EDSTATS" | tee -a $LOG

#Generate total B-factors (the 0-cycle PDB file already has total B-factors)
if ($DOTLS == 0) then
  cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final_tot.pdb
else
  tlsanl \
  XYZIN $WORKDIR/${PDBID}_final.pdb \
  XYZOUT $WORKDIR/${PDBID}_final_tot.pdb \
<<eof > $WORKDIR/mapval.log
    BINPUT t
    BRESID t
    ISOOUT FULL
    NUMERICAL
    END
eof
  if ($status) then
     #TLSANL failed
     echo " o Problem with TLSANL" | tee -a $LOG
     echo "COMMENT: TLSANL cannot calculate total B-factors" >> $DEBUG
     echo "PDB-REDO,$PDBID"                                  >> $DEBUG

     #Use the residual B-factors (not pretty)
     cp $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_final_tot.pdb
  endif
endif

#Run edstats on the original model
edstats.pl -hklin $WORKDIR/${PDBID}_0cyc.mtz -xyzin $WORKDIR/${PDBID}_0cyc.pdb  > $WORKDIR/0cyc_edstats.log

#Reformat the edstats.out file
printf "%s\t%s\t%s\t%s\t%s\n" "RESIDUE" "RSR" "SRSR" "RSCCS" "NGRID" > $WORKDIR/${PDBID}_0cyc.eds
#Use the 'all' or the main chain columns
if (`grep -c 'NPa' $WORKDIR/edstats.out` == 0) then
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "Rm") rsr = i; if ($i == "SRGm") srsr = i; if ($i == "CCSm") rsccs = i; if ($i == "NPm") ngrid = i;}} \
       NR > 1 {printf "%s\t%5.3f\t%5.3f\t%4.2f\t%i\n", $1 "_" $2 "_" $3, $rsr, $srsr, $rsccs, $ngrid}' \
  $WORKDIR/edstats.out >> $WORKDIR/${PDBID}_0cyc.eds
else
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "Ra") rsr = i; if ($i == "SRGa") srsr = i; if ($i == "CCSa") rsccs = i; if ($i == "NPa") ngrid = i;}} \
       NR > 1 {printf "%s\t%5.3f\t%5.3f\t%4.2f\t%i\n", $1 "_" $2 "_" $3, $rsr, $srsr, $rsccs, $ngrid}' \
  $WORKDIR/edstats.out >> $WORKDIR/${PDBID}_0cyc.eds
endif


#Run edstats on the final model
edstats.pl -hklin $WORKDIR/${PDBID}_final.mtz -xyzin $WORKDIR/${PDBID}_final_tot.pdb  > $WORKDIR/final_edstats.log

#Reformat the edstats.out file
printf "%s\t%s\t%s\t%s\t%s\n" "RESIDUE" "RSR" "SRSR" "RSCCS" "NGRID" > $WORKDIR/${PDBID}_final.eds
#Use the 'all' or the main chain columns
if (`grep -c 'NPa' $WORKDIR/edstats.out` == 0) then
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "Rm") rsr = i; if ($i == "SRGm") srsr = i; if ($i == "CCSm") rsccs = i; if ($i == "NPm") ngrid = i;}} \
       NR > 1 {printf "%s\t%5.3f\t%5.3f\t%4.2f\t%i\n", $1 "_" $2 "_" $3, $rsr, $srsr, $rsccs, $ngrid}' \
  $WORKDIR/edstats.out >> $WORKDIR/${PDBID}_final.eds
else
  awk 'NR == 1 {for (i = 1; i <= NF; i++) {if ($i == "Ra") rsr = i; if ($i == "SRGa") srsr = i; if ($i == "CCSa") rsccs = i; if ($i == "NPa") ngrid = i;}} \
       NR > 1 {printf "%s\t%5.3f\t%5.3f\t%4.2f\t%i\n", $1 "_" $2 "_" $3, $rsr, $srsr, $rsccs, $ngrid}' \
  $WORKDIR/edstats.out >> $WORKDIR/${PDBID}_final.eds
endif

#Make plots with R
if ($?RSTAT) then
  $RSTAT --slave --no-save \
  --args $WORKDIR $WORKDIR/${PDBID}_0cyc.eds $WORKDIR/${PDBID}_final.eds $WORKDIR/${PDBID}_dRSCC.svg \
  < $TOOLS/dRSCC.R >>& $WORKDIR/mapval.log
  if ($status) then
    #Error making the figures
    echo " o Problem problem making real-space validation plots" | tee -a $LOG
    echo "COMMENT: error in real-space validation" >> $DEBUG
    echo "PDB-REDO,$PDBID"                         >> $DEBUG
  else
    set PLOTS = 1

    #Get the significant changes
    set RSRB  = `grep 'Better RSR'  $WORKDIR/mapval.log | awk '{print $5}'`
    set RSRW  = `grep 'Worse RSR'   $WORKDIR/mapval.log | awk '{print $5}'`
    set RSCCB = `grep 'Better RSCC' $WORKDIR/mapval.log | awk '{print $4}'`
    set RSCCW = `grep 'Worse RSCC'  $WORKDIR/mapval.log | awk '{print $5}'`

    #Give the summary
    echo " " | tee -a $LOG
    echo '                     Better Worse'  | tee -a $LOG
    echo "Real-space CC      : $RSCCB $RSCCW" | tee -a $LOG
    echo "Real-space R-factor: $RSRB $RSRW"   | tee -a $LOG
    echo " " | tee -a $LOG

    #Compress the plot
    gzip -f -S z $WORKDIR/${PDBID}_dRSCC.svg
  endif
endif

################################################    Ligand validation    #################################################

#Make sure there are ligands
if (`echo $LIG_LIST | grep -c ':'` != 0) then

  echo "****** Ligand validation ******" | tee -a $LOG

  #Create YASARA validation files in parallel
  if ($?YASARA) then
    #Loop over ligands
    foreach LIG (`echo $LIG_LIST | sed 's/:/ /g'`)

      #Specify the residue
      set RESNUM = `echo $LIG | cut -c 2-`
      set CHID   = `echo $LIG | cut -c 1`

      #Return label for job launching
ligvalrunning:

      #Only launch new jobs when the number of cores is not exceeded
      #Strangely direct line counting on the jobs output doesn't work, so a temporary file is needed
      jobs > $WORKDIR/jobs.log

      if (`cat $WORKDIR/jobs.log | wc -l` < $NPROC) then

        echo "-Validating residue $LIG" | tee -a $LOG

        $YASARA -txt $TOOLS/ligval.mcr \
        "oripdb='$WORKDIR/${PDBID}_0cyc.pdb'" \
        "newpdb='$WORKDIR/${PDBID}_final.pdb'" \
        "resnum='$RESNUM'" \
        "chid='$CHID'" > $WORKDIR/ligval_$LIG.log
      else
        #Wait a bit to start again
        sleep 2
        goto ligvalrunning
      endif
    end

    #Wait for the jobs to finish
    wait
  endif


  #Loop over ligands
  foreach LIG (`echo $LIG_LIST | sed 's/:/ /g'`)

    #Specify the residue
    set RESNUM = `echo $LIG | cut -c 2-`
    set CHID   = `echo $LIG | cut -c 1`
    set PDBLIG = `echo "$CHID $RESNUM" | awk '{printf "%s%4d\n", $1, $2}' | sed 's/ /_/g'`
    set RESID  = `grep '^[AH][TE][OT]' $WORKDIR/${PDBID}_final.pdb | cut -c 18-27 | sed 's/ /_/g' | grep $PDBLIG | head -n 1 | cut -c 1-3 | sed 's/_/ /g'`

    #Get real-space values
    set LIGRSRO = `grep "${CHID}_${RESNUM}[[:space:]]" ${PDBID}_0cyc.eds | awk '{print $2}'`
    set LIGCCO  = `grep "${CHID}_${RESNUM}[[:space:]]" ${PDBID}_0cyc.eds | awk '{print $4}'`
    set LIGRSRF = `grep "${CHID}_${RESNUM}[[:space:]]" ${PDBID}_final.eds | awk '{print $2}'`
    set LIGCCF  = `grep "${CHID}_${RESNUM}[[:space:]]" ${PDBID}_final.eds | awk '{print $4}'`

    #Extract YASARA results
    if (-e $WORKDIR/ligval_$LIG.log) then

      #Print atoms with large shifts
      set NSHIFT = `grep 'RMSD over 1 matched atoms' $WORKDIR/ligval_$LIG.log | cut -c 5- | awk -v shiftco=$SHIFTCO '$14 > shiftco {print $5, $4, $2}'| tee $WORKDIR/bigshift_$LIG.log | wc -l`

      #Get RMSD
      set RMSD = `grep '^Residue' $WORKDIR/ligval_$LIG.log | grep RMSD | awk '{print $11}'`

      #Get formation energy
      if (`grep -c 'Energy of formation in object 1' $WORKDIR/ligval_$LIG.log` == 0) then
        set  EFORMO = 'NA'
      else
        set EFORMO = `grep 'Energy of formation in object 1' $WORKDIR/ligval_$LIG.log | cut -d '=' -f 2 | awk '{print $1}'`
      endif
      if (`grep -c 'Energy of formation in object 2' $WORKDIR/ligval_$LIG.log` == 0) then
        set  EFORMF = 'NA'
      else
        set EFORMF = `grep 'Energy of formation in object 2' $WORKDIR/ligval_$LIG.log | cut -d '=' -f 2 | awk '{print $1}'`
      endif

      #Get hydrogen bond energy
      set EHBOO = `grep -A 100 "Start hbond ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hbond ori" | grep Total | awk '{sum += -1*$17} END {print 0.0+sum}'`
      set EHBOF = `grep -A 100 "Start hbond new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hbond new" | grep Total | awk '{sum += -1*$17} END {print 0.0+sum}'`

      #Get number of hydrogen bonds
      set NHBOO = `grep -A 100 "Start hbond ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hbond ori" | grep -c '^Atom'`
      set NHBOF = `grep -A 100 "Start hbond new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hbond new" | grep -c '^Atom'`

      #Get number of bumps
      set NBUMPO = `grep -A 100 "Start bumps ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End bumps ori" | grep Contacts | awk '{sum += $8} END {print 0+sum}'`
      set NBUMPF = `grep -A 100 "Start bumps new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End bumps new" | grep Contacts | awk '{sum += $8} END {print 0+sum}'`

      #Get hydrophobic contacts and strength
      set NHYPHO = `grep -A 100 "Start hydpho ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hydpho ori" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set NHYPHF = `grep -A 100 "Start hydpho new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hydpho new" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set SHYPHO = `grep -A 100 "Start hydpho ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hydpho ori" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`
      set SHYPHF = `grep -A 100 "Start hydpho new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End hydpho new" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`

      #Get pi-pi contacts and strength
      set NPIPIO = `grep -A 100 "Start pipi ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End pipi ori" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set NPIPIF = `grep -A 100 "Start pipi new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End pipi new" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set SPIPIO = `grep -A 100 "Start pipi ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End pipi ori" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`
      set SPIPIF = `grep -A 100 "Start pipi new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End pipi new" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`

      #Get cation-pi contacts and strength
      set NCATPO = `grep -A 100 "Start catpi ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End catpi ori" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set NCATPF = `grep -A 100 "Start catpi new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End catpi new" | grep strength | awk '{sum += $7} END {print 0+sum}'`
      set SCATPO = `grep -A 100 "Start catpi ori" $WORKDIR/ligval_$LIG.log | grep -B 100 "End catpi ori" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`
      set SCATPF = `grep -A 100 "Start catpi new" $WORKDIR/ligval_$LIG.log | grep -B 100 "End catpi new" | grep strength | awk '{sum += $11} END {printf "%5.3f\n", 0+sum}'`
    endif

    #Give summary
    echo " "                                                              | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    echo "****** Ligand validation details ($RESID $CHID $RESNUM) ******" | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    echo " "                                                              | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    echo "                                    Before  Final"       | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    echo "Real-space R-factor               : $LIGRSRO   $LIGRSRF" | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    echo "Real-space correlation            : $LIGCCO    $LIGCCF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    if ($?YASARA) then
      echo "Energy of formation (kJ/mol)      : $EFORMO $EFORMF"   | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Hydrogen bond energy (kJ/mol)     : $EHBOO  $EHBOF"    | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Number of hydrogen bonds          : $NHBOO  $NHBOF"    | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Number of bumps                   : $NBUMPO  $NBUMPF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Number of hydrophobic interactions: $NHYPHO  $NHYPHF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Hydrophobic interaction strength  : $SHYPHO  $SHYPHF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Number of Pi-Pi interactions      : $NPIPIO  $NPIPIF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Pi-Pi interaction strength        : $SPIPIO  $SPIPIF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Number of cation-Pi interactions  : $NCATPO  $NCATPF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Cation-Pi interaction strength    : $SCATPO  $SCATPF"  | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo " "                                                     | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "Atoms shifted more than ${SHIFTCO}A     : $NSHIFT"     | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
      echo "All atom RMSD for residue (A)     : ${RMSD}"           | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
    endif
    echo " "                                                   | tee -a $LOG | tee -a $WORKDIR/${PDBID}_ligval.txt
  end
endif

##############################################   Create quality boxplots  ################################################

#Create a datafile
echo "RFREE,RFFIN,OZRAMA,FZRAMA,OCHI12,FCHI12,URESO"         > $WORKDIR/stats.csv
if ($GOTR == 0 || $ZCALERR == 1) then
  #Compensate for new R-free set
  echo "$RFCALUNB,$RFFIN,$OZRAMA,$FZRAMA,$OCHI12,$FCHI12,$URESO" >> $WORKDIR/stats.csv
else
  echo "$RFCAL,$RFFIN,$OZRAMA,$FZRAMA,$OCHI12,$FCHI12,$URESO" >> $WORKDIR/stats.csv
endif

#Make boxplots with R
if ($?RSTAT) then
  $RSTAT --slave --no-save \
  --args $WORKDIR $TOOLS/pdb_redo_stats.csv $WORKDIR/stats.csv $WORKDIR/${PDBID}_boxplot.svg \
  < $TOOLS/boxplot.R >>& $WORKDIR/boxplot.log
  if ($status) then
    #Error making the figures
    echo " o Problem problem making quality boxplots" | tee -a $LOG
    echo "COMMENT: error in quality boxplots" >> $DEBUG
    echo "PDB-REDO,$PDBID"                    >> $DEBUG
  else
    set PLOTS = 1

    #Compress the plot
    gzip -f -S z $WORKDIR/${PDBID}_boxplot.svg
  endif
endif


################################################   Create YASARA scenes  #################################################

echo " " | tee -a $LOG
echo "****** Creating final output ******" | tee -a $LOG
#Make TLS scene
if ($?YASARA) then
  echo "-Creating YASARA scenes" | tee -a $LOG
  #Make TLS scene
  if ($DOTLS == 1) then
    echo " o TLS-group scene" | tee -a $LOG
    $TOOLS/tls2mcr $WORKDIR/${PDBID}_final.pdb $WORKDIR/${PDBID}_tls.mcr $WORKDIR/${PDBID}_tls.sce > $WORKDIR/yasara.log
    $YASARA -txt $WORKDIR/${PDBID}_tls.mcr >> $WORKDIR/yasara.log
  endif
  
  #Make H-bond scene
  if (-e $WORKDIR/hbond.mcr) then
    echo " o Hydrogen bond restraint scene" | tee -a $LOG
    $YASARA -txt $WORKDIR/hbond.mcr >> $WORKDIR/yasara.log
    mv $WORKDIR/distel_distance_restraints.sce $WORKDIR/${PDBID}_hbond.sce
  endif
  
  #Make KRAB scene
  if (-e $WORKDIR/krab.mcr) then
    echo " o Antibody-specific restraint scene" | tee -a $LOG
    $YASARA -txt $WORKDIR/krab.mcr >> $WORKDIR/yasara.log
    mv $WORKDIR/distel_torsion_restraints.sce $WORKDIR/${PDBID}_torsion.sce
  endif 

  #Make atomic shift scene
  echo " o Atomic shift scene" | tee -a $LOG
  $YASARA -txt $TOOLS/atom_shift.mcr \
  "oripdb='$WORKDIR/${PDBID}_0cyc.pdb'" \
  "newpdb='$WORKDIR/${PDBID}_final.pdb'" \
  "outsce='$WORKDIR/${PDBID}_shift.sce'" >> $WORKDIR/yasara.log

  #Compress the scenes
  tar -cf ${PDBID}_scenes.tar *.sce
  bzip2 -f $WORKDIR/${PDBID}_scenes.tar
else
  echo "-YASARA not installed, cannot make visualisation scenes" | tee -a $LOG
endif

#############################################  Create Refmac command file  ###############################################

#Not for the data bank
if ($LOCAL == 1 || $SERVER == 1) then
  echo "-Creating Refmac command script" | tee -a $LOG

  #Copy all optimised settings
  echo "#Refmac command script from PDB_REDO $VERSION" > $WORKDIR/$PDBID.refmac
  echo "#"                                            >> $WORKDIR/$PDBID.refmac
  echo "#Use of riding hydrogens"                     >> $WORKDIR/$PDBID.refmac
  echo "make hydrogen $HYDROGEN"                      >> $WORKDIR/$PDBID.refmac
  echo "#B-factor model selection"                    >> $WORKDIR/$PDBID.refmac
  echo "refi bref $BREFTYPE"                          >> $WORKDIR/$PDBID.refmac
  if (`echo $REFIRES | grep -c REFI` != 0) then
    echo "#Resolution cutoffs"                        >> $WORKDIR/$PDBID.refmac
    echo "$REFIRES"                                   >> $WORKDIR/$PDBID.refmac
  endif
  echo "#Solvent related settings"                    >> $WORKDIR/$PDBID.refmac
  echo "scal type $SOLVENT $SCALING"                  >> $WORKDIR/$PDBID.refmac
  echo "solvent YES"                                  >> $WORKDIR/$PDBID.refmac
  echo "$MASKPAR"                                     >> $WORKDIR/$PDBID.refmac
  echo "tlsd waters exclude"                          >> $WORKDIR/$PDBID.refmac
  echo "#Restraint weights"                           >> $WORKDIR/$PDBID.refmac
  #Keep the final refinement settings unless there was no final refinement
  if ($RESOTYPE > 3) then
    if ($TLSBEST == "none") then
      set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
    else
      set CWEIGHT = `echo "$WGTSIG MATRIX $TLSBEST"`
    endif
  else if ($BLTBEST == 'auto') then
    set CWEIGHT = `echo "$WGTSIG AUTO 2.50"`
  else
    set CWEIGHT = `echo "$WGTSIG $WGTTYPE $BLTBEST"`
  endif
  echo "weight  $CWEIGHT"                             >> $WORKDIR/$PDBID.refmac
  #Only give the B-factor restraint weight if optimised
  if ($BBEST == "overall" || $BBEST == "none") then
    #Do nothing
  else
    echo "temp $BBEST"                                >> $WORKDIR/$PDBID.refmac
  endif
  if ($TWIN == "twin") then
    echo "#Twinning"                                  >> $WORKDIR/$PDBID.refmac
    echo "$TWIN"                                      >> $WORKDIR/$PDBID.refmac
  endif
  if (`echo $NCSTYPE | grep -c ncs` != 0) then
    echo "#NCS handling"                              >> $WORKDIR/$PDBID.refmac
    echo "$NCSTYPE"                                   >> $WORKDIR/$PDBID.refmac
    echo "$NCSALIGN"                                  >> $WORKDIR/$PDBID.refmac
    echo "$NCSNEIGH"                                  >> $WORKDIR/$PDBID.refmac
    echo "$NCSSTRICT"                                 >> $WORKDIR/$PDBID.refmac
  endif
  if (`echo $JELLY | grep -c ridg` != 0) then
    echo "#Other restraints"                          >> $WORKDIR/$PDBID.refmac
    echo "$JELLY"                                     >> $WORKDIR/$PDBID.refmac
    echo "$TORSION"                                   >> $WORKDIR/$PDBID.refmac
  endif
endif

###################################################   Clean up round 2   #################################################

#Delete FoldX files (silently)
if ($?FOLDX) then
  rm $WORKDIR/rotabase.txt >& /dev/null
  rm $WORKDIR/*.fxout      >& /dev/null
endif

#Create a new version file if it is too old
set NOW  = `date  +%s`
if (-e $TOOLS/pdb_redo.ver) then
  if ($ISOSX == 0) then
    set VERS = `date -r $TOOLS/pdb_redo.ver +%s`
  else
    set VERS = `stat -t %s -f %Sm $TOOLS/pdb_redo.ver`
  endif
else
  set VERS = 0
endif
@ AGE = ($NOW - $VERS)
#Only do this for the databank runs (versions will not change for the server or for local runs)
if ($AGE > 43200 && $LOCAL == 0) then
  #Update reusable versions file or bypass the file creation
  touch $TOOLS/pdb_redo.ver >& /dev/null
  if ($status) then
    #Cannot create reusable versions file
    echo "PDB_REDO     $VERSION" > $WORKDIR/versions.txt
    $TOOLS/versions.csh $TOOLS >> $WORKDIR/versions.txt
  else
    echo "-Updating versions file" | tee -a $LOG
    $TOOLS/versions.csh $TOOLS > $TOOLS/pdb_redo.ver

    #Copy the version data
    echo "PDB_REDO     $VERSION" > $WORKDIR/versions.txt
    cat  $TOOLS/pdb_redo.ver >> $WORKDIR/versions.txt
  endif
else
  #Copy the version data
  echo "PDB_REDO     $VERSION" > $WORKDIR/versions.txt
  cat  $TOOLS/pdb_redo.ver    >> $WORKDIR/versions.txt
endif

#Clean out the mtz files (only for the databank)
if ($LOCAL == 0) then

  #Delete more data?
  set DELLABEL =
  if (`echo $ANOMCOEF | cut -c 1` != "") then
    set DELLABEL = "DELFAN PHDELAN"
  endif

  #Loop over all mtz files
  foreach STAGE (0cyc besttls final)

    #Make a copy
    cp $WORKDIR/${PDBID}_$STAGE.mtz $WORKDIR/allcolumn.mtz

    #Remove the FC PHIC FC_ALL_LS and PHIC_ALL_LS columns
    mtzutils \
    HKLIN  $WORKDIR/allcolumn.mtz \
    HKLOUT $WORKDIR/${PDBID}_$STAGE.mtz \
<<eof >> $WORKDIR/mtzcleanup.log
      EXCLUDE FC PHIC FC_ALL_LS PHIC_ALL_LS $DELLABEL
      END
eof
  end
endif

#################################################    Create webpage    ###################################################
echo "-Consolidating results" | tee -a $LOG

#Calculate percentage
if ($?FRAC) then
  set FRACPC  = `echo $FRAC  | awk '{printf "%4.1f%%", 100*$1}'`
endif

#Create directories and copy files
gzip -f $WORKDIR/${PDBID}_0cyc.pdb
gzip -f $WORKDIR/${PDBID}_0cyc.mtz
gzip -f $WORKDIR/${PDBID}_0cyc.eds
gzip -f $WORKDIR/${PDBID}_besttls.pdb
gzip -f $WORKDIR/${PDBID}_besttls.mtz

set TOUTPUT = "$OUTPUT.tmp"
mkdir -p $TOUTPUT
cp $WORKDIR/${PDBID}_0cyc.pdb.gz          $TOUTPUT/
cp $WORKDIR/${PDBID}_0cyc.mtz.gz          $TOUTPUT/
cp $WORKDIR/${PDBID}_0cyc.eds.gz          $TOUTPUT/
cp $WORKDIR/${PDBID}_besttls.pdb.gz       $TOUTPUT/
cp $WORKDIR/${PDBID}_besttls.mtz.gz       $TOUTPUT/
cp $WORKDIR/${PDBID}_final.pdb            $TOUTPUT/
cp $WORKDIR/${PDBID}_final_tot.pdb        $TOUTPUT/
cp $WORKDIR/${PDBID}_final.mtz            $TOUTPUT/
cp $WORKDIR/${PDBID}_final.eds            $TOUTPUT/
cp $WORKDIR/${PDBID}_final.scm            $TOUTPUT/
cp $WORKDIR/${PDBID}_final.py             $TOUTPUT/
cp $WORKDIR/versions.txt                  $TOUTPUT/
if (-e $WORKDIR/${PDBID}_scenes.tar.bz2) then
  cp $WORKDIR/${PDBID}_scenes.tar.bz2     $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}_ligval.txt) then
  cp $WORKDIR/${PDBID}_ligval.txt         $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}.rtest) then
  cp $WORKDIR/${PDBID}.rtest              $TOUTPUT/
endif
if (-e $WORKDIR/${PDBID}_het.cif) then
  cp $WORKDIR/${PDBID}_het.cif            $TOUTPUT/
  if ($LOCAL == 0) then
    mkdir -p $RESTOUT
    cp $WORKDIR/${PDBID}_het.cif          $RESTOUT/
  endif
endif
if (-e $WORKDIR/renumber.json) then
  cp $WORKDIR/renumber.json               $TOUTPUT/
endif
if ($PLOTS == 1) then
  cp $WORKDIR/${PDBID}_dRSCC.svgz         $TOUTPUT/
  cp $WORKDIR/${PDBID}_boxplot.svgz       $TOUTPUT/
endif
if ($WCERR == 0) then
  mkdir -p $TOUTPUT/wo
  mkdir -p $TOUTPUT/wc
  mkdir -p $TOUTPUT/wf
  cp $WORKDIR/wo/* $TOUTPUT/wo/
  cp $WORKDIR/wc/* $TOUTPUT/wc/
  cp $WORKDIR/wf/* $TOUTPUT/wf/
endif


######################################### Create data files for PDBe sliders ##############################################

#Calculate the quality change scores
set DIV = 3

#Fit to the data
if ($SIGRFFIN == 'NA') then
    set SIGRFU = $SIGRFCAL
else
    set SIGRFU = $SIGRFFIN
endif
if ($TSTCNT < 1 || $GOTR == 0 || $ZCALERR == 1) then
    set ZDFREE = `echo "$RFCALUNB $RFFIN $SIGRFU" | awk '{print ($1 - $2)/$3}'`
else
    set ZDFREE = `echo "$RFCAL $RFFIN $SIGRFU" | awk '{print ($1 - $2)/$3}'`
endif

#Geometric quality
if ($OZPAK2 == 'NA' || $FZPAK2 == 'NA') then
    set DZPAK2 = 0.000
    @ DIV = ($DIV - 1)
else
    set DZPAK2 = `echo "$OZPAK2 $FZPAK2" | awk '{print ($2 - $1)}'`
endif
if ($OZRAMA == 'NA' || $FZRAMA == 'NA') then
    set DZRAMA = 0.000
    @ DIV = ($DIV - 1)
else
    set DZRAMA = `echo "$OZRAMA $FZRAMA" | awk '{print ($2 - $1)}'`
endif
if ($OCHI12 == 'NA' || $FCHI12 == 'NA') then
    set DCHI12 = 0.000
    @ DIV = ($DIV - 1)
else
    set DCHI12 = `echo "$OCHI12 $FCHI12" | awk '{print ($2 - $1)}'`
endif

if ($DIV == 0) then
    set DZSCORE = 'null'
else
    set DZSCORE = `echo "$DZPAK2 $DZRAMA $DCHI12 $DIV" | awk '{print (($1 + $2 + $3)/$4)}'`
endif


#Print the values to the JSON file
echo '{'                                 > $WORKDIR/pdbe.json
echo "  "\""pdbid"\"": "\""$PDBID"\""," >> $WORKDIR/pdbe.json
echo '  "ddatafit": {'                  >> $WORKDIR/pdbe.json
echo "    "\""zdfree"\"": $ZDFREE,"     >> $WORKDIR/pdbe.json
echo '    "range-lower": -12.9,'        >> $WORKDIR/pdbe.json
echo '    "range-upper": 12.9'          >> $WORKDIR/pdbe.json
echo '    },'                           >> $WORKDIR/pdbe.json
echo '  "geometry": {'                  >> $WORKDIR/pdbe.json
echo "    "\""dzscore"\"": $DZSCORE,"   >> $WORKDIR/pdbe.json
echo '    "range-lower": -1.17,'        >> $WORKDIR/pdbe.json
echo '    "range-upper": 1.17'          >> $WORKDIR/pdbe.json
echo '    }'                            >> $WORKDIR/pdbe.json
echo '}'                                >> $WORKDIR/pdbe.json

#Copy over the file
cp $WORKDIR/pdbe.json $TOUTPUT/





########################################### Create summaries for the server ##############################################
#Do the actual writing (note that this is not xhtml, so the code is less formal)
if ($SERVER == 1) then

  #Set values for the summary

  #R-factors as percentages
  if ($RHEAD == 1) then
    set SRFACT = `echo $RFACT | awk '{printf "%4.1f%%", 100*$1}'`
  else
    set SRFACT = "<em>none</em>"
  endif
  if ($RFHEAD == 1) then
    set SRFREE = `echo $RFREE | awk '{printf "%4.1f%%", 100*$1}'`
    set SRFCAL = `echo $RFCAL | awk '{printf "%4.1f%%", 100*$1}'`
  else
    #If there is no R-free given use the 'unbiased' R-free for reference
    set SRFREE = "<em>none</em>"
  endif
  set SRCAL  = `echo $RCAL  | awk '{printf "%4.1f%%", 100*$1}'`
  #Compensate if there is a new R-free set or R-free was severely biased
  if ($GOTR == 0 || $ZCALERR == 1) then
    set SRFCAL = `echo $RFCALUNB | awk '{printf "%4.1f%%", 100*$1}'`
  else
    set SRFCAL = `echo $RFCAL | awk '{printf "%4.1f%%", 100*$1}'`
  endif
  set SRTLS  = `echo $RTLS  | awk '{printf "%4.1f%%", 100*$1}'`
  set SRFTLS = `echo $RFTLS | awk '{printf "%4.1f%%", 100*$1}'`
  set SRFIN  = `echo $RFIN  | awk '{printf "%4.1f%%", 100*$1}'`
  set SRFFIN = `echo $RFFIN | awk '{printf "%4.1f%%", 100*$1}'`
  if ($BBEST == none) then
    set SBBEST = '1.00'
  else
    set SBBEST = $BBEST
  endif

  #Set colours for model changes
  if ($GOTR == 0 || $ZCALERR == 1) then
    set CRFFIN  = `echo $RFFIN $RFCALUNB $SIGRFCAL | awk '{if ($1 < ($2 - 2.6*$3)) {print "better"} else if ($1 > ($2 + 2.6*$3)) {print "worse"} else {print "none"}}'`
  else
    set CRFFIN  = `echo $RFFIN $RFCAL $SIGRFCAL | awk '{if ($1 < ($2 - 2.6*$3)) {print "better"} else if ($1 > ($2 + 2.6*$3)) {print "worse"} else {print "none"}}'`
  endif
  set CCCFFIN = `echo $ZCCF | awk '{if ($1 > 2.60) {print "better"} else if ($1 < -2.60) {print "worse"} else {print "none"}}'`
  set CBRMSZ  = `echo $OBRMSZ $FBRMSZ | awk '{if ($1 > 1.000) {if ($2 < $1) {print "better"} else if ($2 > $1) {print "worse"} else {print "same"}} else if ($1 < 1.000 && $2 > 1.000) {print "worse"} else {print "none"}}'`
  set CARMSZ  = `echo $OARMSZ $FARMSZ | awk '{if ($1 > 1.000) {if ($2 < $1) {print "better"} else if ($2 > $1) {print "worse"} else {print "same"}} else if ($1 < 1.000 && $2 > 1.000) {print "worse"} else {print "none"}}'`

  if ($FZRAMA == 'NA' || $OZRAMA == 'NA') then
    set CZRAMA = 'none'
  else if ($TFZRAMA > $TOZRAMA || $TFZRAMA == 100) then
    set CZRAMA = 'better'
  else if ($TFZRAMA < $TOZRAMA) then
    set CZRAMA = 'worse'
  else
    set CZRAMA = 'none'
  endif
  if ($FCHI12 == 'NA' || $OCHI12 == 'NA') then
    set CCHI12 = 'none'
  else if ($TFCHI12 > $TOCHI12 || $TFCHI12 == 100) then
    set CCHI12 = 'better'
  else if ($TFCHI12 < $TOCHI12) then
    set CCHI12 = 'worse'
  else
    set CCHI12 = 'none'
  endif
  if ($FWBMPS == 'NA' || $OWBMPS == 'NA') then
    set CWBMPS = 'none'
  else if ($TFWBMPS > $TOWBMPS || $TFWBMPS == 100) then
    set CWBMPS = 'better'
  else if ($TFWBMPS < $TOWBMPS) then
    set CWBMPS = 'worse'
  else
    set CWBMPS = 'none'
  endif
  if ($FZPAK2 == 'NA' || $OZPAK2 == 'NA') then
    set CZPAK2 = 'none'
  else if ($TFZPAK2 > $TOZPAK2 || $TFZPAK2 == 100) then
    set CZPAK2 = 'better'
  else if ($TFZPAK2 < $TOZPAK2) then
    set CZPAK2 = 'worse'
  else
    set CZPAK2 = 'none'
  endif
  if ($FHBSAT == 'NA' || $OHBSAT == 'NA') then
    set CHBSAT = 'none'
  else if ($TFHBSAT > $TOHBSAT || $TFHBSAT == 100) then
    set CHBSAT = 'better'
  else if ($TFHBSAT < $TOHBSAT) then
    set CHBSAT = 'worse'
  else
    set CHBSAT = 'none'
  endif

  #Write the summary.log
  #Introduction
  echo "<p>PDB_REDO $VERSION has refined the model in <em>$PDBIN:t</em> against the experimental data in "    > $WORKDIR/summary.log
  if ($HKLIN != "") then
    echo "<em>$HKLIN:t</em> " >> $WORKDIR/summary.log
  else
    echo "<em>$MTZIN:t</em> " >> $WORKDIR/summary.log
  endif
  echo "(which extends to $DATARESH&Aring; resolution and has $REFCNT reflections, $TSTCNT ($TSTPRC%) of "     >> $WORKDIR/summary.log
  echo "which were flagged as a test set). "                                                   >> $WORKDIR/summary.log
  if ($GOTR == 0) then
    echo "A new R-free set had to be created with $FRACPC of the data. "                          >> $WORKDIR/summary.log
  endif
  if ($INREST != "") then
    echo "Geometric restraints in $INREST:t were used during refinement. "                          >> $WORKDIR/summary.log
  endif
  echo "The R and R-free values from the PDB header were $SRFACT and $SRFREE. The recalculated "   >> $WORKDIR/summary.log
  echo "R and R-free values of the starting model were $SRCAL and $SRFCAL.</p>"   >> $WORKDIR/summary.log

  #Twinning and NCS
  if ( ($TWIN == "twin") && ($DONCS == 1)) then
    echo "<p>PDB_REDO detected non-crystallographic symmetry and twinning. "                    >> $WORKDIR/summary.log
  else if ( ($TWIN == "twin") && ($STRICTNCS == 1)) then
    echo "<p>PDB_REDO detected twinning and strict non-crystallographic symmetry. "             >> $WORKDIR/summary.log
  else if ( ($TWIN == "twin") && ($DONCS == 0)) then
    echo "<p>PDB_REDO did not detect non-crystallographic symmetry, but did detect twinning. "  >> $WORKDIR/summary.log
  else if ( ($TWIN != "twin") && ($DONCS == 1)) then
    echo "<p>PDB_REDO detected non-crystallographic symmetry, but did not detect twinning. "    >> $WORKDIR/summary.log
  else if ( ($TWIN != "twin") && ($STRICTNCS == 1)) then
    echo "<p>PDB_REDO detected strict non-crystallographic symmetry, but did not detect twinning. " >> $WORKDIR/summary.log
  else
    echo "<p>PDB_REDO detected neither non-crystallographic symmetry nor twinning. "            >> $WORKDIR/summary.log
  endif

  #Resolution cut-off
  if ($RESOCHECK == 0) then
    echo "The input model was refined using only data to $RESO&Aring; resolution. " >> $WORKDIR/summary.log
    echo '<a href="http://www.cmbi.ru.nl/pdb_redo/faq.html#pair" target="_blank">Paired refinement</a> ' >> $WORKDIR/summary.log
    echo "was used to establish a new resolution cut-off at $URESO&Aring;. $NREFCNT reflections, $NTSTCNT of " >> $WORKDIR/summary.log
    echo "which were flagged as a test set, were used in the rest of the run." >> $WORKDIR/summary.log
  endif

  #B-factors and TLS
  if ($STRICTNCS == 1 && $PMTRIX > 9) then
    echo "As there were $REFPATM reflections per atom and there was strict non-crystallographic symmetry with $PMTRIX copies," >> $WORKDIR/summary.log
  else
    echo "As there were $REFPATM reflections per atom available,"                                  >> $WORKDIR/summary.log
  endif
  echo "$BMODELPHRASE"                                                                           >> $WORKDIR/summary.log
  if ($DOTLS == 1) then
    echo "A TLS model for grouped atom movement "                                    >> $WORKDIR/summary.log
    if ($OPTTLSG == $PDBID) then
      echo "with the TLS groups from input model " >> $WORKDIR/summary.log
    else if ($OPTTLSG == 'REDO') then
      echo "with one TLS group per chain " >> $WORKDIR/summary.log
    else
      echo "with user-defined TLS groups " >> $WORKDIR/summary.log
    endif
    echo "was used. " >> $WORKDIR/summary.log
  endif


  #Score table
  echo "<br><center><h2>Validation metrics before and after PDB_REDO</h2><table>" >> $WORKDIR/summary.log
  echo "<tr><th></th><th>Starting model</th><th>PDB_REDO</th></tr>" >> $WORKDIR/summary.log  #Awful html
  echo "<tr><td>R</td><td class="\""rj"\"">$SRCAL</td><td class="\""rj"\"">$SRFIN</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>R-free</td><td class="\""rj"\"">$SRFCAL</td><td class="\""rj $CRFFIN"\"">$SRFFIN</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Free correlation coefficient &#91;<a href="\""http://www.cmbi.ru.nl/pdb_redo/faq.html#ccfr"\"" target="\""_blank"\"">?</a>&#93;</td><td class="\""rj"\"">$CCFOLD</td><td class="\""rj $CCCFFIN"\"">$CCFFIN</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>rmsZ bonds &#91;<a href="\""http://www.cmbi.ru.nl/pdb_redo/faq.html#rmsz"\"" target="\""_blank"\"">?</a>&#93;</td><td class="\""rj"\"">$OBRMSZ</td><td class="\""rj $CBRMSZ"\"">$FBRMSZ</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>rmsZ angles &#91;<a href="\""http://www.cmbi.ru.nl/pdb_redo/faq.html#rmsz"\"" target="\""_blank"\"">?</a>&#93;</td><td class="\""rj"\"">$OARMSZ</td><td class="\""rj $CARMSZ"\"">$FARMSZ</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Gibbs folding energy (kcal/mol) &#91;<a href="\""http://www.cmbi.ru.nl/pdb_redo/faq.html#fold"\"" target="\""_blank"\"">?</a>&#93;</td><td class="\""rj"\"">$OGFOLD</td><td class="\""rj"\"">$FGFOLD</td></tr>" >> $WORKDIR/summary.log
  echo '<tr><th colspan="3"><center>Percentiles</center></th><tr>' >> $WORKDIR/summary.log
  echo "<tr><td>Ramachandran plot</td><td class="\""rj"\"">$TOZRAMA</td><td class="\""rj $CZRAMA"\"">$TFZRAMA</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Rotamer quality</td><td class="\""rj"\"">$TOCHI12</td><td class="\""rj $CCHI12"\"">$TFCHI12</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Bump severity &#91;<a href="\""http://www.cmbi.ru.nl/pdb_redo/faq.html#wbmp"\"" target="\""_blank"\"">?</a>&#93;</td><td class="\""rj"\"">$TOWBMPS</td><td class="\""rj $CWBMPS"\"">$TFWBMPS</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Fine packing</td><td class="\""rj"\"">$TOZPAK2</td><td class="\""rj $CZPAK2"\"">$TFZPAK2</td></tr>" >> $WORKDIR/summary.log
  echo "<tr><td>Hydrogen bond satisfaction</td><td class="\""rj"\"">$TOHBSAT</td><td class="\""rj $CHBSAT"\"">$TFHBSAT</td></tr></table></center><br>" >> $WORKDIR/summary.log

  #Change table
  echo "<center><h2>Significant model changes</h2><table>" >> $WORKDIR/summary2.log
  echo "<tr><td>Rotamers changed</td><td class="\""rj"\"">$NDROTA</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Side chains flipped</td><td class="\""rj"\"">$HBFLIP</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Side chains completed</td><td class="\""rj"\"">$NSCBLT</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Waters removed</td><td class="\""rj"\"">$NWATDEL</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Peptides flipped</td><td class="\""rj"\"">$NBBFLIP</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Chiralities fixed</td><td class="\""rj"\"">$NCHIRFX</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Residues fitting density better</td><td class="\""rj"\"">$RSCCB</td></tr>" >> $WORKDIR/summary2.log
  echo "<tr><td>Residues fitting density worse</td><td class="\""rj"\"">$RSCCW</td></tr></table></center>" >> $WORKDIR/summary2.log

endif

############################################## Final administrative steps ################################################
#Create summary statistics
echo "$PDBID $VERSION $RFACT $RFREE $RCAL $RFCAL $SIGRFCAL $RFCALUNB $RFCALZ $RTLS $RFTLS $SIGRFTLS $RFTLSUNB $RFTLSZ $RFIN $RFFIN $SIGRFFIN $RFFINUNB $RFFINZ $RFRRAT $BBEST $TLSBEST $BLTBEST $NWATDEL $NBBFLIP $NSCBLT $NDROTA $HBFLIP $STFLIP $NCHIRFX $OZPAK1 $NZPAK1 $FZPAK1 $OZPAK2 $NZPAK2 $FZPAK2 $OZRAMA $NZRAMA $FZRAMA $OCHI12 $NCHI12 $FCHI12 $OBCONF $NBCONF $FBCONF $OBRMSZ $NBRMSZ $FBRMSZ $OARMSZ $NARMSZ $FARMSZ $OBUMPS $NBUMPS $FBUMPS $OHBUNS $NHBUNS $FHBUNS $OGFOLD $NGFOLD $FGFOLD $PROG $DYEAR $RESOLUTION $DATARESH $DATARESL $NREFCNT $TSTCNT $TSTPRC $NTSTCNT $REFPATM $AAXIS $BAXIS $CAXIS $ALPHA $BETA $GAMMA $BAVER $BWILS $BREFTYPE $SOLVENT $VDWPROBE $IONPROBE $RSHRINK $DOTLS $NTLS $OPTTLSG $ORITLS $LEGACY '$SPACEGROUP' $RSCCB $RSCCW $RSRB $RSRW $OWBMPS $NWBMPS $FWBMPS $OHBSAT $NHBSAT $FHBSAT $URESO $CCWOLD $CCWFIN $ZCCW $CCFOLD $CCFFIN $ZCCF $WAVELENGTH $ISTWIN $SOLVD $EXPTYP $COMPLETED $NOPDB $NOSF $USIGMA $ZCALERR $TIME $RESOTYPE $FALSETWIN $TOZRAMA $TFZRAMA $TOCHI12 $TFCHI12 $TOZPAK2 $TFZPAK2 $TOWBMPS $TFWBMPS $TOHBSAT $TFHBSAT $OHRMSZ $NHRMSZ $FHRMSZ $TOZPAK1 $TFZPAK1" > $WORKDIR/data.txt
cp $WORKDIR/data.txt $TOUTPUT/
if ($?COMMENT) then
  $TOOLS/txt2json.py -i $WORKDIR/data.txt -o $TOUTPUT/data.json -c "$COMMENT"
else
  $TOOLS/txt2json.py -i $WORKDIR/data.txt -o $TOUTPUT/data.json
endif

#Give final results (usefull for batch jobs and quick result interpretation)
if ($SERVER == 1) then
  cat $WORKDIR/data.txt >> $BASE/success.txt
endif

#Finish the logfiles and copy them over
if ($LOCAL == 1 || $SERVER == 1) then
  if ($SERVER == 1) then
    echo "</pre>" >> $LOG
    cp $WORKDIR/summary.log  $TOUTPUT/
    cp $WORKDIR/summary2.log $TOUTPUT/
  endif
  #Clean and copy the log file
  grep -v '^\[' $LOG > $TOUTPUT/$PDBID.log
endif

#Copy files over for running jobs in CCP4i or deposition
if ($LOCAL == 1 || $SERVER == 1) then
  cp $WORKDIR/$PDBID.refmac      $TOUTPUT/
  cp $WORKDIR/${PDBID}_final.log $TOUTPUT/
  if ($DOTLS == 1) then
    cp $WORKDIR/${PDBID}_final.tls $TOUTPUT/
  endif
endif


#Swap out directories
if (-e $OUTPUT) then
  mv $OUTPUT $OUTPUT.old 
  mv $TOUTPUT $OUTPUT
  if ($SERVER == 0 && $LOCAL == 1) then
    #Do not remove old directory
  else
    rm -r $OUTPUT.old
  endif
else
  mv $TOUTPUT $OUTPUT
endif

#Finish things for the server
if ($SERVER == 1) then
  #Compress the lot and add a link
  cd $STDIR/output
  zip -r alldata.zip *

  #Delete all unneeded files
  rm ${PDBID}_0cyc.pdb.gz
  rm ${PDBID}_0cyc.mtz.gz
  rm ${PDBID}_besttls.pdb.gz
  rm ${PDBID}_besttls.mtz.gz
  rm ${PDBID}_final.eds
  rm ${PDBID}_final_tot.pdb
  rm ${PDBID}_final.log
  rm versions.txt
  rm ${PDBID}_scenes.tar.bz2
  if (-e ${PDBID}_ligval.txt) then
    rm ${PDBID}_ligval.txt
  endif
  if (-e ${PDBID}.rtest) then
    rm ${PDBID}.rtest
  endif
  if (-e ${PDBID}_het.cif) then
    rm ${PDBID}_het.cif
  endif
  if ($WCERR == 0) then
    rm -rf wo
    rm -rf wc
    rm -rf wf
  endif
  mv data.txt  $STDIR/data.txt

  #Write out status file
  touch $STDIR/processEnded.txt
endif

#Remove tempdir
if ($CLEAN == 1) then
  rm -rf $WORKDIR
endif

echo "-Done. That was a lot of work!" | tee -a $LOG
