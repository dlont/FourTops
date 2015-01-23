import xml.etree.cElementTree as ET
import subprocess
import time
import os
tree = ET.ElementTree(file='config/Run2_Samples.xml')

root = tree.getroot()
datasets = root.find('datasets')
args = []
jobSize = 1000000
for d in datasets:
    if d.attrib['add'] == '1':
        args.append(["./MACRO", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff'], d.attrib['filenames']])
outfiles = []
processes = []
if not os.path.exists("Terminal_Output"):
    os.makedirs("Terminal_Output")

for row in args:
    if row[3] == '1':
        title = row[1]
        execCommands = list(row)
        execCommands.extend(["", ""])
        totalEvents = float(row[8])*float(row[9])
        if (totalEvents > jobSize):
            endEvent = 0
            while (totalEvents-(endEvent*jobSize) > 0):
                startStr = str(endEvent*jobSize)
                endStr = str((endEvent+1)*jobSize)
                execCommands[12] = startStr
                execCommands[13] = endStr
                execCommands[1] = title+"_"+str(endEvent+1)
                outfile = open("Terminal_Output/"+execCommands[1]+".out", 'w')
                outfiles.append(outfile)
                popen = subprocess.Popen(execCommands, stdout=outfile)
                processes.append(popen)
                endEvent += 1
        else:
            execCommands[12] = "0"
            execCommands[13] = str(jobSize)
            outfile = open("Terminal_Output/"+execCommands[1]+".out", 'w')
            outfiles.append(outfile)
            popen = subprocess.Popen(execCommands, stdout=outfile)
            processes.append(popen)
#        popen.wait()
#        for i in row:
#            print i

procsDone = 0
while procsDone < len(processes):
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
    print '{} jobs of {} completed.  Timestamp: {}'.format(procsDone, len(processes), time.ctime())

