from ROOT import TChain
from glob import glob

path = '/pnfs/iihe/cms/store/user/heilman/FOURTOP_NLO_13TeV_aMCatNLO_LHE_pythia8/crab_TOPTREE_742_FOURTOP_1M_v3/150605_093838/0000/*.root'
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
