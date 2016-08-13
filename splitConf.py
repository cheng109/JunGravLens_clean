import numpy as np
import pprint
import subprocess


def readOriginalConf(confName):
    models = {"PTMASS","SIE", "NFW", "SPEMD", "SERSIC" }
    header = ""
    mCombs = []
    f = open(confName, 'r')
    for line in f.readlines():
        items = line.split()
        if len(items)>=1:
            if items[0] not in models and line[0]!='#':
                header += line
            if items[0] in models:
                mCombs.append(line)

    f.close()
    return header, mCombs


def expand(paraFrom, paraTo, paraInc, quota) :
    ret = []


    n = int(np.floor((paraTo - paraFrom)/paraInc))


    if n > quota or quota==0 or n==0:
        ret.append((paraFrom, paraTo, paraInc))
    else:
        for i in range(n):
            ret.append((paraFrom +i*paraInc, paraFrom + (i+1)*paraInc, paraInc))
        #ret.append((paraFrom+(nAvail)*paraInc, paraTo, paraInc))

    quota = int(quota/(n+1))


    return quota, ret


def singleMix(para):
    # for SIE model:  6 levels depth for 6 parameters of SIE model;
    singlePara = []
    for p0  in para[0]:
        for p1 in para[1]:
            for p2 in para[2]:
                for p3 in para[3]:
                    for p4 in para[4]:
                        for p5 in para[5]:
                            singlePara.append([p0[0], p0[1], p0[2],
                                          p1[0], p1[1], p1[2],
                                          p2[0], p2[1], p2[2],
                                          p3[0], p3[1], p3[2],
                                          p4[0], p4[1], p4[2],
                                          p5[0], p5[1], p5[2]])
    return singlePara


def finalMix(mix):
    # depth = 3 ;  for 3 models;
    fMix = []
    for p0 in mix[0]:
        for p1 in mix[1]:
            for p2 in mix[2]:
                fMix.append((p0, p1,p2))
    return fMix

def convertStr(fMix):
    print "Number of final splits: ",  len(fMix)
    fMixConf = []
    for models in fMix:
        oneComb = ""

        for singleM in models:
            tempStr = "SIE "
            for num in singleM:
                tempStr += str(num) + " "
            tempStr += "\n"
            oneComb += tempStr
        fMixConf.append(oneComb)

    return fMixConf

def writeConf(path, header, fMixConf):
    i=0
    fileNameList = []
    for conf in fMixConf:
        fileName =  "testConf_" + str(i)+".txt"
        fileNameList.append(fileName)
        f = open(path + fileName, 'w')
        f.write(header + "\n" + conf )
        f.close()
        i+=1
    return fileNameList

def splitCombs (mCombs):
    mix = []
    # assume all of them are SIE models;
    nLimit = 40
    for model in mCombs:   ## 3 models
        items = model.split()
        i = 1
        para = []
        while i<len(items):
            nLimit, paraList = expand(float(items[i]), float(items[i+1]), float(items[i+2]), nLimit)

            para.append(paraList)
            i += 3
        singlePara = singleMix(para)
        mix.append(singlePara)
    fMix = finalMix(mix)

    fMixConf = convertStr(fMix)
    return fMixConf


def  createRun(runFileName,confFileNameList):
    ind = 0
    f = open(runFileName,'w')
    for fileName in confFileNameList:
        outputFileName = "output_" + fileName
        ind += 1
        if ind%20 == 0 :
            f.write("./junGL horseshoe_test/ " + fileName + " " + outputFileName + " \n"  )
        else:
            f.write("./junGL horseshoe_test/ " + fileName + " " + outputFileName + " &\n"  )
    f.close()
    subprocess.call("chmod +x " + runFileName, shell=True)

def main():

    confName="horseshoe_test/conf.txt"

    header, mCombs = readOriginalConf(confName)
    fMixConf = splitCombs(mCombs)
    confFileNameList = writeConf("horseshoe_test/", header, fMixConf)

    createRun("combRun", confFileNameList)


if __name__=="__main__": 
    main()