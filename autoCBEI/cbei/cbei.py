import sys,re,json,os,shutil
from Bio import SeqIO
from Bio.SeqUtils import GC

baseT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "U": "A",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "u": "a",
    "N": "N",
    "n": "n",
    "R": "Y",
    "Y": "R",
    "M": "K",
    "K": "M",
    "W": "W",
    "S": "S",
    "H": "D",
    "D": "H",
    "B": "V",
    "V": "B",
    "r": "y",
    "y": "r",
    "m": "k",
    "k": "m",
    "w": "w",
    "s": "s",
    "h": "d",
    "d": "h",
    "b": "v",
    "v": "b"
}


def getPAMReg(pam):
    spArr = list(pam)
    for s in range(len(spArr)):
        if (re.match(r"N", spArr[s], re.I)):
            spArr[s] = r"\w"
        elif (re.match(r"R", spArr[s], re.I)):
            spArr[s] = "[AGR]"
        elif (re.match(r"Y", spArr[s], re.I)):
            spArr[s] = "[CTY]"
        elif (re.match(r"M", spArr[s], re.I)):
            spArr[s] = "[ACM]"
        elif (re.match(r"K", spArr[s], re.I)):
            spArr[s] = "[GTK]"
        elif (re.match(r"S", spArr[s], re.I)):
            spArr[s] = "[GCS]"
        elif (re.match(r"W", spArr[s], re.I)):
            spArr[s] = "[ATW]"
        elif (re.match(r"H", spArr[s], re.I)):
            spArr[s] = "[ATCWMYH]"
        elif (re.match(r"B", spArr[s], re.I)):
            spArr[s] = "[GTCKSYB]"
        elif (re.match(r"V", spArr[s], re.I)):
            spArr[s] = "[GACRSMV]"
        elif (re.match(r"D", spArr[s], re.I)):
            spArr[s] = "[GATRKWD]"
        treg = '(?=('+"".join(spArr)+'))'
    return re.compile(treg, re.I)

def tranStrand(seq):
    rSeq = seq[::-1]
    newSeq = ""
    for i in range(len(rSeq)):
        newSeq += baseT[rSeq[i]]
    return newSeq


def changeCodon(codon, pos, change):
    if (pos == 1):
        return change+codon[pos:]
    elif (pos == len(codon)):
        return codon[0:-1]+change
    else:
        return codon[0:pos-1]+change+codon[pos:]


def addLabel(seq, addPos, addStr, reg):
    count = -1
    newSeq = ''
    errLab = 1
    for i in range(len(seq)):
        if (re.match(reg, seq[i])):
            newSeq = newSeq+seq[i]
            continue
        else:
            count += 1
        if (count == addPos):
            newSeq = newSeq+addStr+seq[i]
            errLab = 0
        else:
            newSeq = newSeq+seq[i]
    if (errLab):
        newSeq = newSeq+addStr
    return newSeq


def calSapcer(seq, regPos, beinfo, endReg, mode):
    preStop = []
    pam=beinfo[0]
    spLen=beinfo[1]
    editWin=beinfo[2:4]
    direct=beinfo[4]
    for tpos in regPos:
        pamInfo = {}
        pamInfo['mode'] = mode
        pamPos = tpos.start()+1
        b = -1
        e = -1
        spb = -1
        spe = -1
        if (direct == 5):
            if (pamPos-spLen > -1):
                spb = pamPos-spLen
                spe = pamPos-1
            else:
                continue
            b = spb+editWin[0]-1
            e = spb+editWin[1]-1
        else:
            if (pamPos+spLen > len(seq)):
                continue
            spb = pamPos+len(pam)
            spe = pamPos+len(pam)+spLen-1
            b = pamPos+editWin[0]+len(pam)-1
            e = pamPos+editWin[1]+len(pam)-1
        pamInfo['pamSeq'] = seq[tpos.start():(tpos.start()+len(pam))]
        if (mode == 'Minus'):
            pamInfo['spacerRange'] = str(
                len(seq)-spb+1)+"-"+str(len(seq)-spe+1)
            pamInfo['editRange'] = str(len(seq)-b+1)+"-"+str(len(seq)-e+1)
            pamInfo['pamRange'] = str(
                len(seq)-pamPos+2-len(pam))+"-"+str(len(seq)-pamPos+1)
        else:
            pamInfo['spacerRange'] = str(spb)+"-"+str(spe)
            pamInfo['editRange'] = str(b)+"-"+str(e)
            pamInfo['pamRange'] = str(pamPos)+"-"+str(pamPos+len(pam)-1)

        pamInfo['pureSpacer'] = seq[spb-1:spe]
        # print(pureSpacer)
        longspb = spb
        longspe = spe

        if (spb % 3 == 0):
            longspb = longspb - 2
        elif (spb % 3 == 2):
            longspb = longspb - 1

        if (spe % 3 == 1):
            longspe = longspe + 2
        elif (spe % 3 == 2):
            longspe = longspe + 1

        longSpacer = seq[longspb-1:longspe]
        tedit = []
        tcodonall = []
        for tsp in range(int(len(longSpacer)/3)):
            tcodon = longSpacer[tsp*3:tsp*3+3]
            if (re.match(endReg, tcodon)):
                tedit.append(longspb+tsp*3)
            tcodonall.append(tcodon)
        tcodonORI = tcodonall.copy()

        editPos = ""
        location = -1.1
        editPatten = ""
        tios = []
        for ttedit in tedit:
            teditCodonIndex = int((ttedit-longspb)/3)
            newCodon=""
            for subI in range(3):
                ttIndex = subI+ttedit
                if (ttIndex >= b and ttIndex <= e):
                    if (re.match(r'C', tcodonORI[teditCodonIndex][subI], re.I)):
                        # if (mode == 'Minus'):
                        #     print(tcodonORI[teditCodonIndex]+"\t"+tcodonORI[teditCodonIndex][subI]+"\t"+str(subI))
                        tios.append(ttIndex)
                        newCodon=newCodon+"(c->t)"
                        # tcodonORI[teditCodonIndex][subI]="(c->t)"
                        # tcodonall[teditCodonIndex] = changeCodon(
                        #     tcodonORI[teditCodonIndex], subI+1, "(c->t)")
                    else:
                        newCodon=newCodon+tcodonORI[teditCodonIndex][subI]
                else:
                    newCodon=newCodon+tcodonORI[teditCodonIndex][subI]
            tcodonall[teditCodonIndex]=newCodon

        if (not tios):
            continue

        teditpatten = []
        for ii in range(len(tios)):
            teditpatten.append(seq[tios[ii]-2:tios[ii]])
        editPatten = ",".join(teditpatten)

        if (mode == 'Minus'):
            location = (len(seq)-tios[0]+1) / len(seq)
        else:
            location = tios[0]/len(seq)

        for ios in range(len(tios)):
            if (mode == 'Minus'):
                tios[ios] = str(len(seq)-tios[ios]+1)
            else:
                tios[ios] = str(tios[ios])

        # editPos=str(len(seq)-tios[1]+1)+","+str(len(seq)-tios[0]+1)
        editPos = ",".join(tios)

        if (not editPos):
            continue

        detailSeq = ",".join(tcodonall)

        skipReg = re.compile(r',|\{|\}|\[|\]|\)|\-|\>|[c,t]')

        detailSeq = addLabel(detailSeq, abs(
            longspb-spb)+editWin[0]-1, '[', skipReg)
        detailSeq = addLabel(detailSeq, abs(
            longspb-spb)+editWin[1], ']', skipReg)

        if (direct == 5):
            detailSeq = addLabel(detailSeq, abs(longspb-spb), '{', skipReg)
            detailSeq = addLabel(detailSeq, len(
                longSpacer)-abs(longspe-spe), '}|', skipReg)
            detailSeq = detailSeq+',' + \
                seq[longspe:longspe+len(pam)-(abs(longspe-spe))]
        else:
            detailSeq = addLabel(detailSeq, abs(longspb-spb), '|{', skipReg)
            detailSeq = addLabel(detailSeq, len(
                longSpacer)-abs(longspe-spe)+2, '}', skipReg)
            detailSeq = seq[longspb -
                            len(pam)+abs(longspb-spb)-1:longspb-1]+","+detailSeq
        # print(detailSeq.upper())
        pamInfo['detailSpacer'] = detailSeq.upper()

        pamInfo['editpos'] = editPos
        pamInfo['location'] = location
        pamInfo['editpatten'] = editPatten
        preStop.append(pamInfo)
    return preStop

def runBatch(bename, beinfo, seqDict,rawPath):
    if ( not os.path.exists(rawPath)):
        os.mkdir(rawPath)
    outFileName=bename+"_"+os.path.basename(sys.argv[1])+".cbei";
    endReg = re.compile(r'cag|caa|cga', re.I)
    endRegM = re.compile(r'cca', re.I)
    pamReg = getPAMReg(beinfo[0])
    out=open(os.path.join(rawPath,outFileName),"w");
    for key in seqDict.keys():
        seq = str(seqDict[key].seq)
        seqr = tranStrand(seq)
        pos = pamReg.finditer(seq)
        posR = pamReg.finditer(seqr)
        plusRes = calSapcer(seq, pos, beinfo, endReg, 'Plus')
        minusRes = calSapcer(seqr, posR, beinfo, endRegM, 'Minus')
        plusRes.extend(minusRes)
        plusRes = sorted(plusRes, key=(lambda x: x['location']))
        for tres in plusRes:
            out.write(key+"\t"+"\t".join([
                str(tres['mode']),
                str(tres['location']),
                str(tres['detailSpacer']),
                str(tres['pureSpacer']),
                str(tres['spacerRange']),
                str(tres['editpos']),
                str(tres['editRange']),
                str(tres['pamSeq']),
                str(tres['pamRange']),
                str(tres['editpatten'])
            ])+"\n")
    print("CBEI calculation of base editor \""+bename+"\" , done!")

if __name__ == '__main__':
    beinfos = {
        "BE":["NGG",20,4,8,5],
        "YE1-BE3":["NGG",20,5,7,5],
        "EE-BE3":["NGG",20,5,6,5],
        "YEE-BE3":["NGG",20,6,6,5],
        "VQR-BE3":["NGAN",20,4,11,5],
        "VRER-BE3":["NGCG",20,3,10,5],
        "SaBE":["NNGRRT",21,3,12,5],
        "Sa(KKH)-BE3":["NNNRRT",21,3,12,5],
        "Cas12a–BE":["TTTV",20,10,12,3],
        "Target-AID":["NGG",20,2,4,5],
        "Target-AID-NG":["NG",20,2,4,5],
        "xBE3":["NG",20,4,8,5],
        "BE-PLUS":["NGG",20,4,14,5]
    }

    seqDict = SeqIO.index(sys.argv[1], "fasta")
    for key in beinfos:
        runBatch(key, beinfos[key], seqDict)

    print("Calculate complete!")

