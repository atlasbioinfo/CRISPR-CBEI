import os,re
from autocbei.cbei.autoCBEI import error

def judgeBEF(file):
    rules=["Please follow the rules, eg: .",
                    "BE,NGG,20,4,8,5",
                    "name, PAM, spacer length, edit beg, edit end, PAM at 5' or 3' of spacer"]
    if (os.path.isfile(file)):
        tbeinfo={}
        with open(file,"r") as f:
            for line in f:
                if (len(line.strip())<2):
                    continue
                tmp=line.strip().split(",")
                if(len(tmp)!=6):
                    error(["An incomplete line has been detected in your BE input file (list bellow), please complete it:",line.strip()])
                if (tmp[0] not in tbeinfo):
                    tbeinfo[tmp[0]]=tmp[1:]
                else:
                    error([
                        "A duplicate base editors appears in the file, please check:",
                        tmp[0]+"\t"+",".join(tbeinfo[tmp[0]]),
                        ",".join(tmp)
                    ])
                    
        for tbe in tbeinfo:
            if (not re.match(r'^[ATCGNRYMKSWHBVD]+$',tbeinfo[tbe][0],re.I)):
                error([
                    "The PAM sequence of this record",
                    "\t"+tbe+","+",".join(tbeinfo[tbe]),
                    "in the Base editor file does not match the requirement.",
                    "Please indicate with ATCGNRYMKSWHBVD."
                ])
            if (not (tbeinfo[tbe][1].isdecimal() and 
                    tbeinfo[tbe][2].isdecimal() and
                    tbeinfo[tbe][3].isdecimal() and 
                    tbeinfo[tbe][4].isdecimal())):
                error([
                    "The parameters of base editor:",
                    "\t"+tbe+","+",".join(tbeinfo[tbe]),
                    "has non-integers in it. ",
                    "\n".join(rules)
                ])
            if (int(tbeinfo[tbe][2])>int(tbeinfo[tbe][3])):
                error([
                    "There's an error in the parameters of the basic editor:",
                    "\t"+tbe+","+",".join(tbeinfo[tbe]),
                    "The begin pos: "+tbeinfo[tbe][2]+" is larger than the end"+tbeinfo[tbe][3],
                    "\n".join(rules)
                ])
            if (int(tbeinfo[tbe][2])>=int(tbeinfo[tbe][1]) or int(tbeinfo[tbe][3])>=int(tbeinfo[tbe][1])):
                error([
                    "There's an error in the parameters of the basic editor:",
                    "\t"+tbe+","+",".join(tbeinfo[tbe]),
                    "The begin pos: "+tbeinfo[tbe][2]+" or the end"+tbeinfo[tbe][3]+" is larger than the spacer length "+ str(tbeinfo[tbe][1]),
                    "\n".join(rules)
                ])
            if (not (tbeinfo[tbe][4]=="5" or tbeinfo[tbe][4]=="3")):
                error([
                    "There's an error in the parameters of the basic editor:",
                    "\t"+tbe+","+",".join(tbeinfo[tbe]),
                    "The direction should be 5 or 3, not "+tbeinfo[tbe][4],
                    "\n".join(rules)
                ])
        return tbeinfo               
    error([
        "The "+file+" is not exists, please check!"
    ]) 
        