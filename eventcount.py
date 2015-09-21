from ROOT import TChain
from glob import glob

path = '/pnfs/iihe/cms/store/user/heilman/TopTree/CMSSW_74X_v4/TT_TuneZ2star_13TeV-powheg-pythia6-tauola/crab_TOPTREE-TT_TuneZ2star_13TeV-powheg-pythia6-tauola-RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v3-MINIAODSIM--CMSSW_74X_v4--MCRUN2_74_V9-All/150720_081210/0000/TOPTREE_*.root'
files = glob(path)
root_files = []
for f in files:
	root_files.append('dcap://maite.iihe.ac.be' + f)
print root_files
chain = TChain('eventTree')
for rf in root_files:
	chain.Add(rf)
print 'added files'
nEntries = chain.GetEntries();
print nEntries
