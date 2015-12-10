import xml.etree.cElementTree as ET
import subprocess
import time
import os
import glob
#from __future__ import print_function
#import shutil

tree = ET.ElementTree(file='config/Run2DiLepton_TOPTREES.xml')
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
jobSize = 2000000
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])
        files = ["./MACRO", d.attrib['name'], d.attrib['title'], d.attrib['add'], d.attrib['color'], d.attrib['ls'], d.attrib['lw'], d.attrib['normf'], d.attrib['EqLumi'], d.attrib['xsection'], d.attrib['PreselEff'], d.attrib['filenames']]
        args.append(files)
outfiles = []
fileNames = []
processes = []
tempList = []
if not os.path.exists("Submit_Scripts"):
    os.makedirs("Submit_Scripts")

#for row in args:
#    print "checking args..."
#    if row[3] == '1':
#        title = row[1]
#
#        totalEvents = float(row[8])*float(row[9])
##        tempList = list(row)
##        tempList.extend(["", ""])
#        if (totalEvents > jobSize):
#            endEvent = 0
#            while (totalEvents-(endEvent*jobSize) > 0):
#                startStr = str(endEvent*jobSize)
#                endStr = str((endEvent+1)*jobSize)
#                tempList = list(row)
#                tempList.extend(["", ""])
#                tempList[len(tempList)-2] = startStr
#                tempList[len(tempList)-1] = endStr
#                tempList[1] = title+"_"+str(endEvent+1)
#                fileNames.append("Terminal_Output/"+tempList[1]+".out")
#                execCommands.append(tempList)
#                print 'Job {} Created'.format(tempList[1])
#                #popen = subprocess.Popen(execCommands, stdout=outfile)
#                #processes.append(popen)
#                endEvent += 1
#        else:
#            tempList = list(row)
#            tempList.extend(["", ""])
#            tempList[len(tempList)-2] = "0"
#            tempList[len(tempList)-1] = str(jobSize)
#            fileNames.append("Terminal_Output/"+tempList[1]+".out")
#            execCommands.append(tempList)
#            #popen = subprocess.Popen(execCommands, stdout=outfile)
#            #processes.append(popen)
##        popen.wait()
##        for i in row:
##            print i
for row in args:
    command = row[0]+" "+row[1]+" "+row[2]+" "+row[3]+" "+row[4]+" "+row[5]+" "+row[6]+" "+row[7]+" "+row[8]+" "+row[9]+" "+row[10]+" "+row[11]
    f = open("Submit_Scripts/Submit"+row[1]+".sh",'w')
    f.write('#PBS -q localgrid\n')
    f.write('#PBS -l walltime=02:00:00\n')
    f.write('source /user/heilman/.bash_profile\n')
    f.write('export X509_USER_PROXY=/localgrid/heilman/proxy\n')
    f.write('cd /user/heilman/CMSSW_7_6_0/src/\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('cp -pr ./TopBrussels $TMPDIR/.\n')
#    f.write('cd TopBrussels/FourTops\n')
    f.write('echo $TMPDIR\n')
    f.write('echo ">> Copying work area"\n')

    f.write('cd $TMPDIR/TopBrussels/FourTops\n')
    f.write('ls -l\n')
    f.write(command+'\n')
    f.write('echo ">> Training Complete"\n')
    f.write('ls -l\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/Craneens_* /user/heilman/batchoutput\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/MACRO_* /user/heilman/batchoutput\n')
    f.write('cp -pr $TMPDIR/TopBrussels/FourTops/HistosPtEta*.root /user/heilman/batchoutput\n')

#    popen = subprocess.Popen(row)
#    print 'Job {} begun'.format(row[2])
#    processes.append(popen)
#    procsStarted += 1
#    print 'Jobs {} of {} started.  Timestamp: {}'.format(procsStarted, len(execCommands), time.ctime())
#    while (procsStarted-procsDone) >= (numCores/2):
#        time.sleep(60)
#        procsDone = 0
#        for proc in processes:
#            if proc.poll() != None:
#                procsDone+= 1
#        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())
#while (procsDone != len(execCommands)):  #This loop controls the status output for the last 4 jobs that are still running when the above for loop terminates
#        time.sleep(60)
#        procsDone = 0
#        for proc in processes:
#            if proc.poll() != None:
#                procsDone+= 1
#        print '{} jobs of {} Finished.  Timestamp: {}'.format(procsDone, len(execCommands), time.ctime())



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

