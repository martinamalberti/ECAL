import os
import sys
import numpy
import glob

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-w","--workdir",dest="workdir",type="string",default="LocalCont",help="Name of the directory for jobs")
parser.add_option("-o","--outputname",dest="outputname",type="string",default="graphs",help="Name of the output file")
parser.add_option("-e","--executable",dest="executable",type="string",default="LocalCorractionEB",help="Name of the executable")
parser.add_option("-q","--queue",dest="queue",type="string",default="1nh",help="Name of the queue on lxbatch")
parser.add_option("-r","--minR9",dest="minR9",type="float",default="0.",help="Min R9 value")
parser.add_option("-R","--maxR9",dest="maxR9",type="float",default="999.",help="Max R9 value")
parser.add_option("-b","--minNBC",dest="minNBC",type="float",default="0.",help="Min number of basic clusters")
parser.add_option("-B","--maxNBC",dest="maxNBC",type="float",default="999.",help="Max number of basic clusters")
parser.add_option("-p","--useLocalPhi",dest="useLocalPhi",default=False,action="store_true",help="Analyze local phi")
parser.add_option("-W","--useW",dest="useW",default=False,action="store_true",help="Analyze W events")
parser.add_option("-Z","--useZ",dest="useZ",default=False,action="store_true",help="Analyze Z events")
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true",help="Do not submit jobs")
parser.add_option("","--checkJobs",dest="checkJobs",action="store_true",default=False,help="Checks job status")

(options,args)=parser.parse_args()

path = os.getcwd()
wdir = path+'/'+options.workdir

print path

if not options.checkJobs:
    os.mkdir(wdir)
	
    subLine = './%s -o %s --minR9 %f --maxR9 %f --minNBC %f --maxNBC %f'%(options.executable,options.outputname, options.minR9, options.maxR9, options.minNBC, options.maxNBC)
    if options.useLocalPhi:
        subLine = subLine + ' --useLocalPhi'
    if options.useW:
        subLine = subLine + ' --useW'
    if options.useZ:
        subLine = subLine + ' --useZ'

    print subLine

    os.system('cp %s/%s %s'%(path,options.executable,wdir))   
    jobscript = open('%s/sub.sh'%(wdir),'w')
    jobscript.write('cd %s \n'%(wdir))
    jobscript.write('export SCRAM_ARCH=slc5_amd64_gcc462 \n')
    jobscript.write('eval ` scramv1 runtime -sh ` \n')
    jobscript.write('if ( \n')
    jobscript.write('\t touch %s/sub.run \n'%(wdir))
    jobscript.write('\t %s \n'%(subLine))
    jobscript.write(') then \n')
    jobscript.write('\t mv ./%s.root %s \n'%(options.outputname,wdir))
    jobscript.write('\t rm %s/sub.run \n'%(wdir))
    jobscript.write('\t touch %s/sub.done \n'%(wdir))
    jobscript.write('else \n')
    jobscript.write('\t touch %s/sub.fail \n'%(wdir))
    jobscript.write('fi \n')
    os.system('chmod a+x %s/sub.sh'%(wdir))
    print 'bsub -q %s -o %s/sub.log %s'%(options.queue,wdir,jobscript.name )
    if not options.dryRun:
        os.system('bsub -q %s -o %s/sub.log %s'%(options.queue,wdir,jobscript.name ))
       
else:
    if options.workdir == '':
        print 'Please specify for which task you want to check the status. Use option: --workdir'
    else:
        if os.path.isfile('%s/sub.run'%(wdir)):
            print 'job in status: RUNNING'
        elif os.path.isfile('%s/sub.done'%(wdir)):    
            print 'job in status: DONE'
        elif os.path.isfile('%s/sub.fail'%(wdir)):    
            print 'job in status: FAILED'
            print 'You can resubmit with:'
            print 'bsub -q %s -o %s/sub.log %s'%(options.queue,wdir,jobscript.name )
