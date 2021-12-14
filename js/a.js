/*
Apache License
Version 2.0, January 2004

Copyright 2019 Yu Haopeng atlasbioin4@gmail.com atlas_hp@163.com

Licensed under the Apache License, Version 2.0 (the "License");

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*

This file is responsible for alignment, that is, sequence alignment using Web worker locally. I know that the current algorithm is not optimal and only work.

I hope we can work together to improve this algorithm and improve the efficiency of the front-end alignment.

------Yours Haopeng Yu
*/
let header = new Object();
let resultInfoSet = new Set();
onmessage = function (e) {
    let fileTemp = e.data[0];
    let spacerInfo = e.data[1];
    let misNum = e.data[2];
    let reg = e.data[3];
    let subSeq = e.data[4];
    if (subSeq) {
        subSeq = fastaTrim(subSeq)
    }
    runAlign(fileTemp, spacerInfo, misNum, reg, subSeq)
};

function runAlign(fileTemp, spacerInfo, misNum, reg, subSeq) {
    let total = fileTemp.size;
    let step = 1024 * 1024 * 50;
    let gobackstep = 100;
    let reader = new FileReader();
    let cuLoaded = 0;
    let goback = 0;
    reader.onload = async function (e) {
        let seq = fastaTrim(reader.result);
        if (seq[seq.length - 1].fasHeader) {
            header.header = seq[seq.length - 1].fasHeader;
            header.len = seq[seq.length - 1].seq.length
        }
        let tmpPer = parseInt((cuLoaded + step) / total * 10000) / 100;
        if (tmpPer > 100) {
            tmpPer = 100
        };
        postMessage([1, "genomePro", tmpPer]);
        postMessage([2, "genomePro", "Selected Fasta file processing " + tmpPer + "%"]);
        let tmpReturn = await runSubAlign(seq, spacerInfo, misNum, reg, tmpPer, goback);
        if (tmpPer == 100) {
            postMessage([1, "genomePro", 100]);
            postMessage([2, "genomePro", "Calculation complete. Preparing for display"]);
            postMessage([3]);
            postMessage(resultInfoSet)
        }
        console.log(tmpReturn);
        cuLoaded += e.loaded;
        if (cuLoaded < total) {
            readBlob(cuLoaded - gobackstep);
            goback = 1
        } else {
            cuLoaded = total
        }
    };
    reader.onerror = function () {
        console.log(reader.error)
    };
    postMessage([1, "genomePro", 0]);
    postMessage([2, "genomePro", "Selected Fasta file processing"]);
    readBlob(0);

    function readBlob(start) {
        let blob = fileTemp.slice(start, start + step);
        reader.readAsText(blob)
    }
}

function runSubAlign(seq, spacerInfo, misNum, reg, tmpPer, goback) {
    return new Promise(function (resolve, reject) {
        for (let sp = 0; sp < spacerInfo.length; sp++) {
            let misInfo = autoAlign(spacerInfo[sp], seq, misNum, reg);
            postMessage([1, "pro", ((sp + 1) * 100 / spacerInfo.length)]);
            postMessage([2, "pro", "Spacer " + spacerInfo[sp].index + " alignment is completed. The next spacer in alignment"]);
            if (sp == spacerInfo.length - 1) {
                postMessage([1, "pro", 100]);
                postMessage([2, "pro", "All spacer complete for this part of genome"])
            }
            if (misInfo.length == 0) {
                resultInfoSet.add(0 + "BOSSA" + spacerInfo[sp].index)
            } else {
                for (let q = 0; q < misInfo.length; q++) {
                    resultInfoSet.add(1 + "BOSSA" + misInfo[q].index + "BOSSA" + misInfo.length + "BOSSA" + spacerInfo[0].updown + "BOSSA" + misInfo[q].mispos + "BOSSA" + misInfo[q].seq + "BOSSA" + misInfo[q].strand + "BOSSA" + misInfo[q].header + "BOSSA" + misInfo[q].pos + "BOSSA" + misInfo[q].newPAM)
                }
            }
        }
        resolve(tmpPer + " OK")
    })
}

function union(setA, setB) {
    let _union = new Set(setA);
    for (let elem of setB) {
        _union.add(elem)
    }
    return _union
}

function fastaTrim(seq) {
    let seqArr = seq.split("\n");
    let fastaInfo = new Array();
    let headerReg = /^>/;
    let ntCode = /[^atcgn]/i;
    let zhushi = /^\#/;
    let tempSeq = "";
    let headInfo = new Object();
    for (let i = 0; i < seqArr.length; i++) {
        let temp = seqArr[i].replace(/[\r|\n]/g, "");
        if (zhushi.test(temp)) {
            continue
        }
        if (headerReg.test(temp)) {
            temp = temp.replace(/>/, "");
            if (tempSeq !== "") {
                if (!headInfo.fasHeader) {
                    headInfo.fasHeader = header.header;
                    headInfo.len = header.len + tempSeq.length;
                    headInfo.seq = tempSeq
                } else {
                    headInfo.seq = tempSeq;
                    headInfo.len = tempSeq.length
                }
                fastaInfo.push(headInfo);
                headInfo = {};
                tempSeq = ""
            }
            headInfo.fasHeader = temp
        } else {
            if (ntCode.test(temp)) {
                header.header = header.header + temp
            } else {
                tempSeq += temp
            }
        }
    }
    if (!headInfo.fasHeader) {
        headInfo.fasHeader = header.header;
        headInfo.seq = tempSeq;
        headInfo.len = header.len + tempSeq.length
    } else {
        headInfo.seq = tempSeq;
        headInfo.len = tempSeq.length
    }
    fastaInfo.push(headInfo);
    headInfo = {};
    seqArr = [];
    return fastaInfo
}

function rmSubSeq(subSeq, seq) {
    let subSeqR = transStrand(subSeq);
    for (let j = 0; j < seq.length; j++) {
        let pos = seq[j].seq.search(subSeq);
        let posr = seq[j].seq.search(subSeqR);
        if (pos > -1) {
            seq[j].seq = (seq[j].seq.substr(0, pos)) + (seq[j].seq.substr(pos + subSeq.length));
            return seq
        } else if (posr > -1) {
            seq[j].seq = (seq[j].seq.substr(0, posr)) + (seq[j].seq.substr(posr + subSeqR.length));
            return seq
        }
    }
    return seq
}

function autoAlign(spacerInfo, seq, misNum, reg) {
    let sp = spacerInfo.spacer;
    let spr = transStrand(spacerInfo.spacer);
    let lenPAM = spacerInfo.PAM.length;
    let updown = spacerInfo.updown;
    let index = spacerInfo.index;
    if (misNum == 0) {
        return ala0(sp, spr, seq, reg, lenPAM, updown, index)
    } else {
        return alaplus(sp, spr, seq, reg, lenPAM, updown, index, misNum)
    }
}

function spSpacer(spacer, num) {
    let a = parseInt(spacer.length / num);
    let subSpacer = new Array();
    for (let i = 0; i < num - 1; i++) {
        let temp = new Object();
        temp.seq = spacer.substr(i * a, a);
        temp.len = i * a;
        subSpacer.push(temp)
    }
    let temp = new Object();
    temp.seq = spacer.substr((num - 1) * a);
    temp.len = (num - 1) * a;
    subSpacer.push(temp);
    return subSpacer
}

function alaplus(sp, spr, seq, reg, lenPAM, updown, index, misNum) {
    let misInfo = new Array();
    let subsp = spSpacer(sp, parseInt(misNum) + 1);
    let sps = sp.toLowerCase().split("");
    let subspr = spSpacer(spr, parseInt(misNum) + 1);
    let spsr = spr.toLowerCase().split("");
    let misSet = new Set();
    for (let j = 0; j < seq.length; j++) {
        for (let i = 0; i < subsp.length; i++) {
            let spReg = new RegExp(subsp[i].seq, "ig");
            while (spReg.exec(seq[j].seq) !== null) {
                let pos = spReg.lastIndex - subsp[i].seq.length - subsp[i].len;
                let newPAM = "";
                if (updown == "+") {
                    newPAM = seq[j].seq.substr(pos + sp.length, lenPAM)
                } else {
                    newPAM = seq[j].seq.substr(pos - lenPAM, lenPAM)
                }
                if (!reg.tpam.test(newPAM) || pos < 0 || (seq[j].seq.length - pos < sp.length + 1) || (seq[j].seq.length < sp.length + 1)) {
                    continue
                }
                let newSeq = al(seq[j].seq.substr(pos, sp.length), sps, misNum, subsp[i].seq.length, subsp[i].len);
                if (newSeq) {
                    if (seq[j].len != seq[j].seq.length) {
                        pos = seq[j].len - seq[j].seq.length + pos + 1
                    } else {
                        pos++
                    }
                    misSet.add(newSeq.misSeq + "BOSSA" + newSeq.posInfo + "BOSSA+BOSSA" + index + "BOSSA" + seq[j].fasHeader + "BOSSA" + pos + "BOSSA" + newPAM)
                }
            }
            let sprReg = new RegExp(subspr[i].seq, "ig");
            while (sprReg.exec(seq[j].seq) !== null) {
                let pos = sprReg.lastIndex - subspr[i].seq.length - subspr[i].len;
                let newPAM = "";
                if (updown == "+") {
                    newPAM = seq[j].seq.substr(pos - lenPAM, lenPAM)
                } else {
                    newPAM = seq[j].seq.substr(pos + spr.length, lenPAM)
                }
                if (!reg.tpamr.test(newPAM) || pos < 0 || (seq[j].seq.length - pos < sp.length + 1) || (seq[j].seq.length < sp.length + 1)) {
                    continue
                }
                let newSeq = al(seq[j].seq.substr(pos, spr.length), spsr, misNum);
                if (newSeq) {
                    if (seq[j].len != seq[j].seq.length) {
                        pos = seq[j].len - seq[j].seq.length + pos + 1
                    } else {
                        pos++
                    }
                    misSet.add(newSeq.misSeq + "BOSSA" + newSeq.posInfo + "BOSSA-BOSSA" + index + "BOSSA" + seq[j].fasHeader + "BOSSA" + pos + "BOSSA" + newPAM)
                }
            }
        }
    }
    for (let tmp of misSet) {
        let info = tmp.split("BOSSA");
        let tmpInfo = new Object();
        tmpInfo.seq = info[0];
        tmpInfo.strand = info[2];
        tmpInfo.index = info[3];
        tmpInfo.mispos = info[1];
        tmpInfo.header = info[4];
        tmpInfo.pos = info[5];
        tmpInfo.newPAM = info[6];
        misInfo.push(tmpInfo)
    }
    return misInfo
}

function ala0(sp, spr, seq, reg, lenPAM, updown, index) {
    let spReg = new RegExp(sp, "ig");
    let sprReg = new RegExp(spr, "ig");
    let misInfo = new Array();
    for (let j = 0; j < seq.length; j++) {
        while (spReg.exec(seq[j].seq) !== null) {
            reg.tpam.test("YUHAOPENG");
            let pos = spReg.lastIndex - sp.length;
            let newPAM = "";
            if (updown == "+") {
                newPAM = seq[j].seq.substr(spReg.lastIndex, lenPAM)
            } else {
                newPAM = seq[j].seq.substr(pos - lenPAM, lenPAM)
            }
            if (!reg.tpam.test(newPAM)) {
                continue
            }
            if (seq[j].len != seq[j].seq.length) {
                pos = seq[j].len - seq[j].seq.length + pos + 1
            } else {
                pos++
            }
            let tmpInfo = new Object();
            tmpInfo.seq = sp;
            tmpInfo.strand = "+";
            tmpInfo.index = index;
            tmpInfo.pos = pos;
            tmpInfo.header = seq[j].fasHeader;
            tmpInfo.newPAM = newPAM;
            misInfo.push(tmpInfo)
        }
        while (sprReg.exec(seq[j].seq) !== null) {
            reg.tpamr.test("YUHAOPENG");
            let pos = sprReg.lastIndex - sp.length;
            let newPAM = "";
            if (updown == "+") {
                newPAM = seq[j].seq.substr(pos - lenPAM, lenPAM)
            } else {
                newPAM = seq[j].seq.substr(spReg.lastIndex, lenPAM)
            }
            if (!reg.tpamr.test(newPAM)) {
                continue
            }
            if (seq[j].len != seq[j].seq.length) {
                pos = seq[j].len - seq[j].seq.length + pos + 1
            } else {
                pos++
            }
            let tmpInfo = new Object();
            tmpInfo.seq = spr;
            tmpInfo.strand = "-";
            tmpInfo.index = index;
            tmpInfo.pos = pos;
            tmpInfo.header = seq[j].fasHeader;
            tmpInfo.newPAM = newPAM;
            misInfo.push(tmpInfo)
        }
    }
    return misInfo
}

function al(tseq, sps, misNum, partLen, len) {
    tseq = tseq.toLowerCase();
    let arr = tseq.split("");
    let mis = 0;
    let nSeq = "";
    let pos = "";
    for (let i = 0; i < arr.length; i++) {
        if (arr[i] !== sps[i]) {
            pos += (i + 1) + ",";
            mis++;
            nSeq += arr[i].toUpperCase()
        } else {
            nSeq += arr[i]
        }
        if (mis > misNum) {
            return 0
        }
    }
    let align = {
        "misSeq": nSeq,
        "posInfo": pos
    };
    return align
}

function termiSudden() {
    close()
}

function transStrand(seq) {
    "use strict";
    let baseT = {
        A: "T",
        T: "A",
        C: "G",
        G: "C",
        U: "A",
        a: "t",
        t: "a",
        c: "g",
        g: "c",
        u: "a",
        ",": ",",
        "-": "-",
        "[": "]",
        "]": "["
    };
    let rSeq = seq.split("").reverse();
    for (let i = 0; i < rSeq.length; i++) {
        rSeq[i] = baseT[rSeq[i]]
    }
    return rSeq.join("")
}