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
./SLMACRObTag TTJets t\bar{t}+jets Madgraph MLM 1 633 1 2 1 13500.9 831.76 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2-CMSSW_74X_v8-MCRUN2_74_V9/151020_160750/0000/TOPTREE_*.root
echo ">> Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/MACRO_* /user/lbeck/batchoutput
