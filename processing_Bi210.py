import os 
import glob

fileindexlist=[]
indexlist=[]

#for file in glob.glob("/data/snoplus/OfficialProcessing/production_5_0/Solar_5.0.1/Bi210/root/SolarBi210_r*"):
for file in glob.glob("/data/snoplus/OfficialProcessing/production_5_3_0/Bi210Full/SolarBi210_r*"):
    start = file.find('_r')+2
    end = file.find('_s')
    print(file[start:end])
    indexlist.append(file[start:end])
    fileindexlist.append(file)


#for x in range(0,len(fileindexlist)):
for x in range(0,20):

	inFile = open('processing_classifiers.mac', 'r')
	outFile = open('macs/processing_classifiers_Bi210_'+str(x)+'.mac', 'w')
	
	inFileStr=inFile.read()

	inroot= "/rat/inroot/load "+fileindexlist[x]+"\n"

	outFile.write(inroot)

	inFileStr=inFileStr.replace("/rat/procset file \"test.ntuple.root\"","/rat/procset file \"SolarBi210."+str(indexlist[x])+".ntuple.root\"")
	outFile.write(inFileStr)

	inFile.close()
	outFile.close()
	inFile=""

	inFile = open('processing_temp.scr', 'r')
	outFile = open('scripts/processing_Bi210_'+str(x)+'.scr', 'w')

	inFileStr=inFile.read()

	inFileStr= inFileStr.replace("rat -l /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/logs/ /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/processing_Bipo212_46.mac","rat -l /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/logs/Bi210_"+str(x)+".log -o Bi210_"+str(x)+".ntuple.root /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/macs/processing_classifiers_Bi210_"+str(x)+".mac") 
	inFileStr=inFileStr.replace("cp BiPo212_FromOP_140.ntuple.root /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/output/","cp Bi210_"+str(x)+".ntuple.root /data/snoplus/liggins/year1/fitting/alphaBetaClassifier/output2/")
	outFile.write(inFileStr)
#
#	
#
#
	inFile.close()
	outFile.close()
	inFileStr=""
#
	os.system("source /opt/sge/default/common/settings.sh")
	os.system("qsub -cwd -q snoplusSL6 scripts/processing_Bi210_"+str(x)+".scr")
#

