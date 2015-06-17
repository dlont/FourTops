from ROOT import TChain
from glob import glob

path = '/pnfs/iihe/cms/store/user/fblekman/TopTree/TopTree74v1/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_TOPTREE-WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1-MINIAODSIM--TopTree74v1--MCRUN2_74_V9-All/150610_100409/0000/*.root'
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
