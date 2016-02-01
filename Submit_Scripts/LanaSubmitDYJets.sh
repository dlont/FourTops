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
./SLMACRO DYJets DY+jets 1 214 1 2 1 801.15 6025.2 0.0 /pnfs/iihe/cms/store/user/fblekman/TopTree/CMSSW_74X_v8/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v3-CMSSW_74X_v8-MCRUN2_74_V9/151029_124746/0000/TOPTREE_*.root
echo ">> Complete"
ls -l
cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput
cp -pr $TMPDIR/TopBrussels/FourTops/SLMACRO_* /user/lbeck/batchoutput
