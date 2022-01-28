import sys,re,json,os,shutil
import argparse

def error(sentences):
    print("-"*30)
    print("Error:")
    print("\n".join(sentences))
    print("autoCBEI.py exit...")
    print("-"*30)
    sys.exit()

def mainAutoCBEI():
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



    try:
        from Bio import SeqIO
    except ImportError:
        error(["The \"biopython\" package was not found.",
            "Please install \"biopython\"."])

    try:
        import matplotlib
    except ImportError:
        error(["The \"matplotlib\" package was not found.",
            "Please install \"matplotlib\"."])

    try:
        from autocbei.cbei import cbei,stat
    except ImportError:
        error(["The \"cbei\" package was not found.",
            "Please reinstall \"autocbei\" with \"pip install --upgrade autocbei\"."])
    
    from autocbei.cbei.judges import judgeBEF
    # demoPath=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),"Bacillus_subtilis.part500.cds.all.fa")

    parser = argparse.ArgumentParser(description='Enter fasta file of CDSs,\
            output base editor\'s potential editing site and statistics information.\n \
            The demo CDS Fasta file: "[PATHON_HOME]/site-packages/autocbei/Bacillus_subtilis.part500.cds.all.fa".')
    parser.add_argument("cds",metavar='CDS.fasta',help="CDSs in fasta format.",type=str)
    parser.add_argument("-bef",help="[Optional] The file containing the base editor definition, seperate with \",\". eg: \n \
                        \"BE,NGG,20,4,8,5\", that is, \"name, PAM, spacer length, edit beg, edit end, PAM at 5' or 3' of spacer\". \
                        13 base editors are included by default. \
                        ",type=str)
    # parser.add_argument("-be",help="[Optional] Define a base editor, seperate with \",\". eg: \n \
                        # \"BE,NGG,20,4,8,5\", that is, \"name, PAM, spacer length, edit beg, edit end, PAM at 5' or 3' of spacer\"\
                        # ",type=str)
                        
    parser.add_argument("-ns","--nostat",help="Only run CBEI design without statistics and plot.",action="store_true")
    parser.add_argument("-o","--outprefix",help="Directory prefixes can be customized. Default: \"CBEI\" (CBEIRaw, CBEIPlot, CBEIRes). ", type=str)
    args = parser.parse_args()
    # 1. The potential editing sites for the base editor are calculated, 
    #     and each base editor generates a separate cbei file.

    pre="CBEI"
    if (args.outprefix):
        pre=args.outprefix
    
    if (args.bef):
        beinfos=judgeBEF(args.bef)
        print("Successfully imported "+str(len(beinfos))+" base editors from file: "+args.bef)


    rawPath=pre+"Raw"
    plotPath=pre+"Plot"
    resPath=pre+"Stat"

    print("Input file: "+args.cds)
    print("Output directory: "+rawPath)
    print("Base editors: ")
    print("\tBE\tPAM\tSpacer\tEditBegin\tEditEnd\tDirection")
    for key in beinfos.keys():
        print("\t"+key,end="\t")
        for t in beinfos[key]:
            print(t,end="\t")
        print()

    seqDict = SeqIO.index(args.cds, "fasta")
    tName=os.path.basename(args.cds).split(".")
    fileName=tName[0]
    for key in beinfos:
        print("Start calculating: "+key)
        cbei.runBatch(key, beinfos[key], seqDict,rawPath)
    print("Calculate complete!")

    if (not args.nostat):
        print("Begin statistics...")
        print("The statistics directory: "+resPath)
        print("The plot directory: "+plotPath)

        stat.statCBEI(seqDict,rawPath,plotPath,resPath,fileName)
        print("CBEI statistics complete")

if __name__ == "__main__":
    mainAutoCBEI()