from ROOT import TChain
from glob import glob

path = '/pnfs/iihe/cms/store/user/heilman/Skimmed-TopTrees/CMSSW_74X_v7/Muon/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext3-v1-CMSSW_74X_v7-MCRUN2_74_V9/151015_224531/0000/TOPTREE_*.root'
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
