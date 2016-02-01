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
./SLMACRO ttttNLO NLO SM tttt 1 40 1 2 0.563529176 80869000 .008823 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1-CMSSW_74X_v8-MCRUN2_74_V9/151020_160817/0000/TOPTREE_*.root
echo ">> Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/SLMACRO_* /user/lbeck/batchoutput
