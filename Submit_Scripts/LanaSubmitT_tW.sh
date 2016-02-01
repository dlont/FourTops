#PBS -q localgrid
#PBS -l walltime=10:00:00
source /user/lbeck/.bash_profile
export X509_USER_PROXY=/localgrid/lbeck/proxy
cd /user/lbeck/CMSSW_7_6_0/src/
eval `scramv1 runtime -sh`
cp -pr ./TopBrussels $TMPDIR/.
cp -pr /user/lbeck/lib/* $TMPDIR/.
cp -pr /user/lbeck/lib/* /localgrid/lbeck/lib
echo $TMPDIR
echo ">> Copying work area"
cd $TMPDIR/TopBrussels/FourTops
ls -l
./SLMACRO T_tW t tW 1 7 1 2 1 10168.8 35.6 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_160546/0000/TOPTREE_*.root
echo ">> Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/SLMACRO_* /user/lbeck/batchoutput
