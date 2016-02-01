#PBS -q localgrid
#PBS -l walltime=01:00:00
source /user/lbeck/.bash_profile
export X509_USER_PROXY=/localgrid/lbeck/proxy
cd /user/lbeck/CMSSW_7_6_0/src/
eval `scramv1 runtime -sh`
cp -pr ./TopBrussels $TMPDIR/.
echo $TMPDIR
echo ">> Copying work area"
cd $TMPDIR/TopBrussels/FourTops
ls -l
./SLMACRO Data_Run2015D_v3_Silver Data 1 1 1 2 1 1 8202260. 0.0 /pnfs/iihe/cms/store/user/lbeck/Skimmed-TopTrees/CMSSW_74X_v8/Muon/SingleMuon/crab_SingleMuon-Run2015D-05Oct2015-v1-CMSSW_74X_v8-SILVER-74X_dataRun2_Prompt_v2_/02122015_193141/Skimmed_TOPTREE_*.root
echo ">> Training Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/MACRO_* /user/lbeck/batchoutput
