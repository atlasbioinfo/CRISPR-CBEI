import sys,re,json,os,shutil
import argparse

#Set the Base editor parameter
# [PAM, spacer length, edit beg, edit end, direction]
# Direction refers to spacer at the 5' or 3' end of PAM sequence (5 or 3, respectively).
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

rawPath="./CBEIRaw"
plotPath="./CBEIPlot"
resPath="./CBEIRes"

try:
    from cbei import cbei,stat
except ImportError:
    print("Error:\n\t The \"cbei\" package was not found.\n\tPlease copy the \"cbei\" folder to the current directory and run again")
    print("autoCBEI.py exit...")
    sys.exit()

try:
    from Bio import SeqIO
except ImportError:
    print("Error:\n\tThe \"biopython\" package was not found.\n\tPlease install \"biopython\" using \"pip install biopython\"")
    print("autoCBEI.py exit...")
    sys.exit()

parser = argparse.ArgumentParser(description='Enter fasta file of transcripts,\
        output base editor\'s potential editing site and statistics information.')
parser.add_argument("trans",help="Transcripts in fasta format",type=str)
args = parser.parse_args()

# 1. The potential editing sites for the base editor are calculated, 
#     and each base editor generates a separate cbei file.

print("Input file: "+sys.argv[1])
print("Base editors: ")
print("\tBE\tPAM\tSpacer\tEditBegin\tEditEnd\tDirection")
for key in beinfos.keys():
    print("\t"+key,end="\t")
    for t in beinfos[key]:
        print(t,end="\t")
    print()

seqDict = SeqIO.index(sys.argv[1], "fasta")
tName=os.path.basename(sys.argv[1]).split(".")
fileName=tName[0]
for key in beinfos:
    print("Start calculating: "+key)
    cbei.runBatch(key, beinfos[key], seqDict,rawPath)
print("Calculate complete!")
print("Begin statistics...")
stat.statCBEI(seqDict,rawPath,plotPath,resPath,fileName)
print("CBEI statistics complete")