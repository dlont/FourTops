import xml.etree.cElementTree as ET
import subprocess
import time
tree = ET.ElementTree(file='config/Run2SingleLepton_samples.xml')
#tree = ET.ElementTree(file='config/Run2_Samples.xml')

root = tree.getroot()
datasets = root.find('datasets')
args = []
for d in datasets:
    if d.attrib['add'] == '1':
        args.append(["./SLMACRO", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff'], d.attrib['filenames']])

outfiles = []
processes = []

for row in args:
    if row[3] == '1':
        outfile = open(row[1]+".out", 'w')
        outfiles.append(outfile)
        popen = subprocess.Popen(row, stdout=outfile)
        processes.append(popen)
#        popen.wait()
#        for i in row:
#            print i

procsDone = 0
while procsDone < len(args):
    time.sleep(60)
    procsDone = 0
    counter = 0
    for proc in processes:
        counter += 1
        if proc.poll() != None:
            procsDone += 1
            print 'Job {} complete.'.format(counter)
        else:
            print 'Job {} still running.'.format(counter)
    print '{} jobs of {} completed.  Timestamp: {}'.format(procsDone, len(args), time.ctime())
    
