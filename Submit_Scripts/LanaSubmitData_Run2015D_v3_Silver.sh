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
./SLMACRO Data_Run2015D_v3_Silver Data 1 1 1 2 1 1 8202260. 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8-SILVER/SingleMuon/crab_SingleMuon-Run2015D-05Oct2015-v1-CMSSW_74X_v8-SILVER-74X_dataRun2_Prompt_v2/151126_112419/0000/TOPTREE_*.root
echo ">> Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/SLMACRO_* /user/lbeck/batchoutput
