import os,re,random
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
from Bio.SeqUtils import GC
from cycler import cycler

def statCBEI(seqDict,cbeiPath,plotPath,statPath,fileName):
    transNum=len(seqDict)
    if (not os.path.exists(plotPath)):
        os.mkdir(plotPath)

    if (not os.path.exists(statPath)):
        os.mkdir(statPath)

    files=os.listdir(cbeiPath)
    allRaio={}
    for cfile in files:
        if(re.match(r'.*'+fileName+'.*\.cbei$',cfile,re.I)):
            tmp=cfile.split(".")
            be=tmp[0].split("_")[0]
            filePath=os.path.join(cbeiPath,cfile)
            cbeiInfo=getInfo(be,fileName,filePath,statPath)
            allRaio[be]=[len(cbeiInfo["t25"])/transNum,len(cbeiInfo["t5"])/transNum,len(cbeiInfo["t75"])/transNum,]        
            editablePie(cfile,transNum,cbeiInfo["t25"],cbeiInfo["t5"],cbeiInfo["t75"],plotPath)
    print("#"*20)
    print("Pie charts for different BEs have been generated. \nPath:")
    print("\t"+os.path.join(plotPath,"[BE names]_"+fileName+".statPie.png"))
    
    transGCStat(fileName,seqDict,plotPath)
    print("Transcript statistics and mapping completed.\nPath:\n\t"+os.path.join(plotPath,fileName+"_CDSlength.png"))
    print("\t"+os.path.join(plotPath,fileName+"_GCstats.png"))
    print("\t"+os.path.join(plotPath,fileName+"_CodonUsage.png"))
    print("\t"+os.path.join(plotPath,fileName+"_CDSlength.png"))

    statAllRatio(allRaio,plotPath,fileName)
    print("The comparison of CBEI ratio of different BE has been completed.\nPath:")
    print("\t"+os.path.join(plotPath,fileName+".statBar.png"))
    begROC(transNum, cbeiPath, plotPath,fileName)
    print("\t"+os.path.join(plotPath,fileName+".statROC.png"))
    
def begROC(transNum, cbeiPath, plotPath,fileName):
    beName={}
    files=os.listdir(cbeiPath)
    for file in files:
        if not (re.match(r".*"+fileName+".*\.cbei$",file,re.I)):
            continue
        title=file.strip().split(".")
        beName[title[0]]=[0 for k in range(0,100)]
        with open(os.path.join(cbeiPath,file)) as f:
            tgene={}
            for line in f:
                if (re.match(r'^#.*$',line,re.I)):
                    continue
                tsp=line.strip().split("\t")
                if (tsp[0] in tgene):
                    if (tsp[2]<tgene[tsp[0]]):
                        tgene[tsp[0]]=tsp[2]
                else:
                    tgene[tsp[0]]=tsp[2]
            for key in tgene.keys():
                tpos=int(float(tgene[key])*100)+1
                for p in range(tpos,100):
                    beName[title[0]][p]+=1
    statROC(beName,plotPath,fileName,transNum)

def to_percent(value, pos):
    return '%1.0f'%(100*value) + '%'

def statROC(beName,plotPath,fileName,transNum):
    nxarr= [[ 0 for j in range(0,100)]for i in beName.keys()]
    nyarr= [[ 0 for j in range(0,100)]for i in beName.keys()]
    thre=[ i/100 for i in range(0,100)]
    index=0
    names=[]
    for key in beName.keys():
        tname=key.split("_")
        names.append(tname[0])
        for i in range(0,100):
            beName[key][i]/=transNum
            nxarr[index][i]=thre[i]
            nyarr[index][i]=beName[key][i]
        index+=1
    plt.close()
    fig,ax=plt.subplots(figsize=(8,6))
    tlines=['-', '--', ':', '-.']
    for i in range(len(nxarr)):
        ax.plot(nxarr[i],nyarr[i],label=names[i],linestyle=tlines[random.randint(0,3)])
    ax.legend()
    ax.set_xlabel("Threshold")
    ax.set_ylabel("CBEI ratio")
    plt.gca().yaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.gca().xaxis.set_major_formatter(FuncFormatter(to_percent))
    plt.savefig(os.path.join(plotPath,fileName+".statROC.png"),dpi=300)
    plt.close()
    

def statAllRatio(ratio,plotPath,fileName):
    thres=['0.25', '0.5', '0.75']
    xlab=[]
    thre25=[]
    thre5=[]
    thre75=[]
    x = np.arange(len(ratio))
    fig2, ax2 = plt.subplots(3,1, figsize=(8,12),sharex=True)

    for key in ratio:
        xlab.append(key)
        thre25.append(ratio[key][0])
        thre5.append(ratio[key][1])
        thre75.append(ratio[key][2])

    ax2[0].bar(xlab, thre25)
    ax2[0].set(title="Threshold=25%")
    ax2[0].set_xticklabels(xlab,rotation=45)
    autolabel(ax2[0],xlab,thre25)
    ax2[0].set_ylabel("CBEI ratio")
    ax2[1].bar(xlab, thre5)
    ax2[1].set_xticklabels(xlab,rotation=45)
    ax2[1].set_ylabel("CBEI ratio")
    autolabel(ax2[1],xlab,thre5)
    ax2[1].set(title="Threshold=50%")
    ax2[2].bar(xlab, thre75)
    autolabel(ax2[2],xlab,thre75)
    ax2[2].set_xticklabels(xlab,rotation=45)
    ax2[2].set_ylabel("CBEI ratio")
    ax2[2].set(title="Threshold=75%")
    ax2[0].get_yaxis().set_major_formatter(FuncFormatter(to_percent))
    ax2[1].get_yaxis().set_major_formatter(FuncFormatter(to_percent))
    ax2[2].get_yaxis().set_major_formatter(FuncFormatter(to_percent))
    plt.tight_layout()
    plt.savefig(os.path.join(plotPath,fileName+".statBar.png"),dpi=300)
    plt.close() 

def autolabel(plt,x,y):
    for xx, yy in zip(x,y):
        plt.text(xx, yy+0.01, "%1.2f" % (yy*100)+"%", ha='center')

def func(pct, allvals):
    absolute = int(pct/100.*np.sum(allvals))
    return "{:.2f}%\n{:d}".format(pct, absolute)


def editablePie(cfile, transNum,t25,t5,t75,plotPath):
    tmp=cfile.split(".")
    fileName="_".join(tmp[0].split("_")[1:])
    be=tmp[0].split("_")[0]
    labels = 'Editable', 'Others'
    explode = (0.1, 0)
    fig, axs = plt.subplots(1,3,figsize=(9, 4))
    # fig.subtitle=cfile
    data25=[len(t25),transNum-len(t25)]
    data5=[len(t5),transNum-len(t5)]
    data75=[len(t75),transNum-len(t75)]
    # dataBase25=statBase(t25)
    # dataBase5=statBase(t5)
    # dataBase75=statBase(t75)

    axs[0].pie(data25, explode=explode, labels=labels,autopct=lambda pct: func(pct, data25),
            shadow=False, startangle=90, textprops={'size': 'smaller'})
    axs[0].set(aspect="equal", title='Threshold=0.25')
    axs[1].pie(data5, explode=explode, labels=labels,autopct=lambda pct: func(pct, data5),
            shadow=False, startangle=90, textprops={'size': 'smaller'})
    axs[1].set(aspect="equal", title='Threshold=0.5')
    axs[2].pie(data75, explode=explode, labels=labels,autopct=lambda pct: func(pct, data75),
            shadow=False, startangle=90, textprops={'size': 'smaller'})
    axs[2].set(aspect="equal", title='Threshold=0.75')
    plt.suptitle("Editable transcripts ratio of "+be+" of "+fileName)
    plt.savefig(os.path.join(plotPath,be+"_"+fileName+".statPie.png"),dpi=300)
    plt.close()
    
def transGCStat(fileName,seqDict,plotPath):
    import re
    Codons = {
        "TTT": 0, "TTC": 0, "TTA": 0, "TTG": 0,
        "CTT": 0, "CTC": 0, "CTA": 0, "CTG": 0,
        "ATT": 0, "ATC": 0, "ATA": 0, "ATG": 0,
        "GTT": 0, "GTC": 0, "GTA": 0, "GTG": 0,
        "TAT": 0, "TAC": 0, "TAA": 0, "TAG": 0,
        "CAT": 0, "CAC": 0, "CAA": 0, "CAG": 0,
        "AAT": 0, "AAC": 0, "AAA": 0, "AAG": 0,
        "GAT": 0, "GAC": 0, "GAA": 0, "GAG": 0,
        "TCT": 0, "TCC": 0, "TCA": 0, "TCG": 0,
        "CCT": 0, "CCC": 0, "CCA": 0, "CCG": 0,
        "ACT": 0, "ACC": 0, "ACA": 0, "ACG": 0,
        "GCT": 0, "GCC": 0, "GCA": 0, "GCG": 0,
        "TGT": 0, "TGC": 0, "TGA": 0, "TGG": 0,
        "CGT": 0, "CGC": 0, "CGA": 0, "CGG": 0,
        "AGT": 0, "AGC": 0, "AGA": 0, "AGG": 0,
        "GGT": 0, "GGC": 0, "GGA": 0, "GGG": 0}
    gc_values=[]
    seqlength=[]
    cCount=0
    for key in seqDict.keys():
        seq = str(seqDict[key].seq)
        if len(seq) % 3 !=0 :
            continue
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            if (codon in Codons):
                Codons[codon]+=1
                cCount+=1
        seqlength.append(len(seq))
        gc_values.append(GC(seq))
    gc_values=sorted(gc_values)
    plt.subplots(figsize=(6, 6))
    plt.plot(gc_values)
    plt.title("GC%% of %i CDSs \n from %.2f%% to %.2f%%" % (len(seqDict),min(gc_values),max(gc_values)))
    plt.ylabel("GC content (%)")
    plt.xlabel("CDSs")
    plt.savefig(os.path.join(plotPath,fileName+"_GCstats.png"),dpi=300)
    plt.close()
    seqlength=sorted(seqlength)
    plt.plot(seqlength)
    plt.title("CDSs length of %i CDSs \n from %int to %int" % (len(seqDict),min(seqlength),max(seqlength)))
    plt.ylabel("Sequene length (nt)")
    plt.xlabel("CDSs")
    plt.savefig(os.path.join(plotPath,fileName+"_CDSlength.png"),dpi=300)
    plt.close()
    codonName=[]
    codonFre=[]
    tcolor=[]
    for key in Codons.keys():
        Codons[key]=Codons[key]/cCount*100
        codonName.append(key)
        codonFre.append(Codons[key])
        if (re.match(r'CAG|CGA|CAA|CCA',key,re.I)):
            tcolor.append("tomato")
        else:
            tcolor.append("dodgerblue")
    fig3, ax3 = plt.subplots(figsize=(20, 5))
    ax3.bar(codonName, codonFre,color=tcolor)
    ax3.set(title="Codon frequency (%%) \n Frequency of CAG,CGA,CAA,CCA correspond to %.2f%%,%.2f%%,%.2f%%,%.2f%%" % (Codons["CAG"],Codons["CGA"],Codons["CAA"],Codons["CCA"]))
    ax3.set_xticklabels(codonName,rotation=45)
    plt.ylabel("Codon frequency (%)")
    plt.savefig(os.path.join(plotPath,fileName+"_CodonUsage.png"),dpi=300)
    plt.close()


def getInfo(be,fileName,filePath,statPath):
    cbeiTran25={}
    cbeiTran5={}
    cbeiTran75={}
    outf25=open(os.path.join(statPath,be+"_"+fileName+".Thre025.tsv"),"w")
    outf5=open(os.path.join(statPath,be+"_"+fileName+".Thre05.tsv"),"w")
    outf75=open(os.path.join(statPath,be+"_"+fileName+".Thre075.tsv"),"w")
    outf25.write("#Editable sites in the first 25% of the transcript\n")
    outf5.write("#Editable sites in the first 50% of the transcript\n")
    outf75.write("#Editable sites in the first 75% of the transcript\n")
    fheader="\t".join([
        "#Transcript",
        "Strand",
        "Position",
        "CBEIdetail",
        "Spacer",
        "SpacerRegion",
        "EditPosition",
        "EditWindowsRegion",
        "PAM",
        "PAMregion",
        "Pattern"
    ])
    outf25.write(fheader+"\n")
    outf5.write(fheader+"\n")
    outf75.write(fheader+"\n")
    with open(filePath,"r") as tcbei:
        for line in tcbei:
            tarr=line.strip().split("\t")
            tarr[2]=float(tarr[2])
            if tarr[2]<0.25:
                if tarr[0] in cbeiTran25:
                    continue
                else:
                    cbeiTran25[tarr[0]]=tarr[10]
                    outf25.write(line)
                if tarr[0] in cbeiTran5:
                    continue
                else:
                    cbeiTran5[tarr[0]]=tarr[10]
                    outf5.write(line)
                if tarr[0] in cbeiTran75:
                    continue
                else:
                    cbeiTran75[tarr[0]]=tarr[10]
                    outf75.write(line)
            elif tarr[2]<0.5:
                if tarr[0] in cbeiTran5:
                    continue
                else:
                    cbeiTran5[tarr[0]]=tarr[10]
                    outf5.write(line)
                if tarr[0] in cbeiTran75:
                    continue
                else:
                    cbeiTran75[tarr[0]]=tarr[10]
                    outf75.write(line)
            elif tarr[2]<0.75:
                if tarr[0] in cbeiTran75:
                    continue
                else:
                    cbeiTran75[tarr[0]]=tarr[10]
                    outf75.write(line)
    return {"t25":cbeiTran25,"t5":cbeiTran5,"t75":cbeiTran75}

if __name__ == "__main__":
    
    from Bio import SeqIO
    seqDict = SeqIO.index("../Bacillus_subtilis.ASM69118v1.cds.all.fa", "fasta")
    statCBEI(seqDict,"../CBEIRaw/","../CBEIPlot/","../CBEIRes")
