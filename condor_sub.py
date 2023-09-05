#!/usr/bin/env python
#import os, re, sys, math, time, calendar
import os, re, sys, commands, math, time, calendar

print '\nSTART\n'

ts = calendar.timegm(time.gmtime())
fileName = "analysisEtaMuMuGamma.root"
jobName = "analysisEtaMuMuGamma_MINIAOD_1"
jobScript = "cmsRun.sh"
jobCfg = "standardData.py"
rel = "CMSSW_13_0_6"
#eosDir = "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/" + os.environ["USER"] + "/condor/" + jobName + "_" + str(ts) + "/"
eosDir = "/afs/cern.ch/work/e/elfontan/private/ScoutingParkingPaper/CMSSW_13_0_6/src/Run3ScoutingAnalysisTools/" + jobName + "_" + str(ts) + "/"
rootDir = os.environ["CMSSW_BASE"] + "/src/Run3ScoutingAnalysisTools/"
jobDir = rootDir + jobName + "_" + str(ts) + "/"
ret = 0

print("rootDir",rootDir)

fileList = rootDir + "bennet_EtaToMuMuGamma_Run3_2022_MINIAOD_1.list" 
print("fileList",fileList)
nEvents = 1400000                                                                    
nJobs = 700     
#nEvents = 7000000    
#nJobs = 3499    

while ret == 0:
   ret = os.system("mkdir " + jobDir)
   ret = os.system("mkdir " + jobDir + "out/")
   ret = os.system("mkdir " + jobDir + "err/")
   ret = os.system("mkdir " + jobDir + "log/")
   ret = os.system("mkdir " + eosDir)
   ret = os.chdir(os.environ["CMSSW_BASE"]+"/../")
   ret = os.system("tar --exclude='analysisEtaMuMuGamma.root' --exclude='ignore' --exclude='.git' " + 
                   "-zcf " + jobName + ".tgz " + rel)
   ret = os.system("mv " + jobName + ".tgz " + eosDir) 
   ret = os.chdir(rootDir)

   with open(jobDir + jobName + '.jdl', 'w') as jdl:
      jdl.write("universe = vanilla\n")
      jdl.write("x509userproxy = $ENV(X509_USER_PROXY)\n")
      jdl.write("Executable = " + jobScript + "\n")
      jdl.write("Should_Transfer_Files = YES\n")
      jdl.write("WhenToTransferOutput = ON_EXIT\n")
      jdl.write("Transfer_Input_Files = " + jobScript + ", " + jobCfg + "\n")
      jdl.write("Output = "    + jobDir + "out/$(ProcId).out\n")
      jdl.write("Error = "     + jobDir + "err/$(ProcId).err\n")
      jdl.write("Log = "       + jobDir + "log/$(ProcId).log\n")
      jdl.write("Arguments = " + eosDir + " " + jobName + " " + rel + " " + 
                " $(ProcId) " + str(nEvents) + " " + jobCfg + " " + fileList + " " + fileName + "\n")
      jdl.write("+MaxRuntime = 28800\n")
      jdl.write("on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n")
      jdl.write("max_retries = 3\n")
      jdl.write("requirements = Machine =!= LastRemoteHost\n")
      jdl.write("Queue " + str(nJobs) + "\n")      

   os.system("condor_submit " + jobDir + jobName + ".jdl")
   print str(nJobs) + " jobs submitted."
   print "\nYour jobs:"
   os.system("condor_q")
   print
   sys.exit(0)

print("Submission failed.")
sys.exit(1)
