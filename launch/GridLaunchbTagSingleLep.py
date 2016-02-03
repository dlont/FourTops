import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
#from __future__ import print_function
#import shutil

tree = ET.ElementTree(file='config/BtagHistoTTJets.xml')
#tree = ET.ElementTree(file='config/Testing.xml')

root = tree.getroot()
datasets = root.find('datasets')

print "found  "  + str(len(datasets)) + " datasets"

procsDone = 0
procsStarted = 0
numCores = 8
args = []
execCommands = []
topTrees = []
jobSize = 200000
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        files = ["./SLMACRObTag", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff'], d.attrib['filenames']]
        args.append(files)
outfiles = []
fileNames = []
processes = []
tempList = []
if not os.path.exists("Submit_Scripts"):
    os.makedirs("Submit_Scripts")

for row in args:
    command = row[0]+" "+row[1]+" "+row[2]+" "+row[3]+" "+row[4]+" "+row[5]+" "+row[6]+" "+row[7]+" "+row[8]+" "+row[9]+" "+row[10]+" "+row[11]
    f = open("Submit_Scripts/Submit"+row[1]+".sh",'w')
    f.write('#PBS -q localgrid\n')
    f.write('#PBS -l walltime=01:00:00\n')
    f.write('source /user/lbeck/.bash_profile\n')
    f.write('export X509_USER_PROXY=/localgrid/lbeck/proxy\n')
    f.write('cd /user/lbeck/CMSSW_7_6_0/src/\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('cp -pr ./TopBrussels $TMPDIR/.\n')
#    f.write('cd TopBrussels/FourTops\n')
    f.write('echo $TMPDIR\n')
    f.write('echo ">> Copying work area"\n')

    f.write('cd $TMPDIR/TopBrussels/FourTops\n')
    f.write('ls -l\n')
    f.write(command+'\n')
    f.write('echo ">> Complete"\n')
    f.write('ls -l\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/lbeck/batchoutput\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/Histo_* /user/lbeck/batchoutput\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/MACRO_* /user/lbeck/batchoutput\n')

