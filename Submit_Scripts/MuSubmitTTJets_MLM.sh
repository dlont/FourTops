#PBS -q localgrid
#PBS -l walltime=09:00:00
source /user/lbeck/.bash_profile
export X509_USER_PROXY=/localgrid/lbeck/proxy
cd /user/lbeck/CMSSW_7_6_0/src/
eval `scramv1 runtime -sh`
cp -pr ./TopBrussels $TMPDIR/.
echo $TMPDIR
echo ">> Copying work area"
cd $TMPDIR/TopBrussels/FourTops
ls -l
./SLMACROLocal TTJets_MLM t\bar{t}+jets_Madgraph_MLM 1 633 1 2 1 5825.7271 831.76 0.0 /pnfs/iihe/cms/store/user/lbeck/Skimmed-TopTrees/CMSSW_74X_v8/Muon/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/02122015_214501/Skimmed_TOPTREE_*.root
echo ">> Training Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/MACRO_* /user/lbeck/batchoutput
