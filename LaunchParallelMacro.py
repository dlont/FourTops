import xml.etree.cElementTree as ET
import subprocess
import time
import os

tree = ET.ElementTree(file='config/Run2DiLepton_TOPTREES.xml')

root = tree.getroot()
datasets = root.find('datasets')
procsDone = 0
procsStarted = 0
numCores = 8
args = []
execCommands = []
jobSize = 1000000
for d in datasets:
    if d.attrib['add'] == '1':
        args.append(["./MACRO", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff'], d.attrib['filenames']])
outfiles = []
fileNames = []
processes = []
tempList = []
if not os.path.exists("Terminal_Output"):
    os.makedirs("Terminal_Output")

for row in args:
    if row[3] == '1':
        title = row[1]

        totalEvents = float(row[8])*float(row[9])
        if (totalEvents > jobSize):
            endEvent = 0
            while (totalEvents-(endEvent*jobSize) > 0):
                startStr = str(endEvent*jobSize)
                endStr = str((endEvent+1)*jobSize)
                tempList = list(row)
                tempList.extend(["", ""])
                tempList[12] = startStr
                tempList[13] = endStr
                tempList[1] = title+"_"+str(endEvent+1)
                fileNames.append("Terminal_Output/"+tempList[1]+".out")
                execCommands.append(tempList)
                print 'Job {} Created'.format(tempList[1])
                #popen = subprocess.Popen(execCommands, stdout=outfile)
                #processes.append(popen)
                endEvent += 1
        else:
            tempList[12] = "0"
            tempList[13] = str(jobSize)
            fileNames.append("Terminal_Output/"+tempList[1]+".out")
            execCommands.append(tempList)
            #popen = subprocess.Popen(execCommands, stdout=outfile)
            #processes.append(popen)
#        popen.wait()
#        for i in row:
#            print i
for i, row in enumerate(execCommands):
    outfile = open(fileNames[i], 'w')
    print 'File {} opened'.format(fileNames[i])
    outfiles.append(outfile)
    popen = subprocess.Popen(row, stdout=outfiles[i])
    print 'Job {} begun'.format(row[1])
    processes.append(popen)
    procsStarted += 1
    print 'Jobs {} of {} started.  Timestamp: {}'.format(procsStarted, len(execCommands), time.ctime())
    while (procsStarted-procsDone) >= (numCores/2):
        time.sleep(60)
        procsDone = 0
        for proc in processes:
            if proc.poll() != None:
                procsDone+= 1
        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())
while (procsDone != len(execCommands)):  #This loop controls the status output for the last 4 jobs that are still running when the above for loop terminates
        time.sleep(60)
        procsDone = 0
        for proc in processes:
            if proc.poll() != None:
                procsDone+= 1
        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())



#while procsDone < len(processes):
#    time.sleep(60)
#    procsDone = 0
#    counter = 0
#    for proc in processes:
#        counter += 1
#        if proc.poll() != None:
#            procsDone += 1
#            print 'Job {} complete.'.format(counter)
#        else:
#            print 'Job {} still running.'.format(counter)
#    print '{} jobs of {} completed.  Timestamp: {}'.format(procsDone, len(processes), time.ctime())

