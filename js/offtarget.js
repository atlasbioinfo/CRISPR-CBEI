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

function judgeFileFirst(fileName) {
    let fileTemp = document.getElementById(fileName).files[0];
    if (fileTemp === undefined) {
        alert("Please select the local fasta file for off-target prediction." + '\n' + "This file can be a genomes sequence or multiple orfs, but in fasta format that is not compressed.");
        return 0;
    } else {
        return fileTemp;
    }
}


function offPredict2() {

    $("#offResultDir").show();
    let fileTemp = judgeFileFirst("fastaUpload2");
    if (!fileTemp) {
        return;
    }
    let begin = new Date();
    let spArr = getSpacerText("spacerSequence2");
    let misNum = $("#misNum2").val();
    let spacerInfo = new Array();
    let infoDraw = getValueForDraw2();

    for (let i = 0; i < spArr.length; i++) {
        if (spArr[i] == "") {
            continue
        };
        let tempInfo = new Object();
        tempInfo.spacer = spArr[i];
        tempInfo.index = i + 1;
        tempInfo.PAM = infoDraw[1];
        tempInfo.updown = infoDraw[6];
        spacerInfo.push(tempInfo);
    }
    let pam = getRegPAM(infoDraw[1]);
    pam.tpam = new RegExp(pam.tpam, 'i');
    pam.tpamr = new RegExp(pam.tpamr, 'i');
    let subSeq = "";
    $("#proBar").attr('style', "width: 5%;");
    $("#genomeProBar").attr('style', "width: 5%;");
    $("#proInfo").text("Data processing, please wait");
    $("#genomeProInfo").text("Reading genome file");
    //这里的subSeq意思的排除所预测的ORF序列，这个需要在网页上再加一个框，或者点击Off-target预测后直接复制过来
    if (typeof (Worker) !== "undefined") {
        myWorker = new Worker('js/align.js');
        $("#myModal").modal("show");
        myWorker.postMessage([fileTemp, spacerInfo, misNum, pam, subSeq]);
        myWorker.onmessage = function (e) {
            if (e.data[0] == 1) {
                $("#" + e.data[1] + "Bar").attr("style", "width: " + e.data[2] + "%;");
            } else if (e.data[0] == 2) {
                $("#" + e.data[1] + "Info").text(e.data[2]);
            } else if (e.data[0] == 3) {
                closeModal();
                let end = new Date();
                let time = end - begin;
                console.log("Time consumption is " + time + "ms");
            } else {
                generaOffTableSet(e.data, spArr, "offListBody2");
            }
        }
    } else {
        showSomething("warnNotSupport");
    }
}


function getSpacerText(textareaName) {
    let upSpacer = $("#" + textareaName).val();
    let spArr = upSpacer.split("\n");
    let newArr = new Array();
    for (let i = 0; i < spArr.length; i++) {
        if (/[^atcg\s\r\n]/i.test(spArr[i])) {
            continue;
        }
        if (/\s|\r|\n/.test(spArr[i])) {
            spArr[i] = spArr[i].replace(/\s|\r|\n/g, "");
        }
        if (spArr[i].length < 16) {
            continue;
        }
        newArr.push(spArr[i]);
    }
    return newArr;
}

// function offPredict() {
//     let fileTemp = judgeFileFirst("fastaUpload");
//     if (!fileTemp) { return; }
//     if (!judgeCheckBox()) {
//         return;
//     }
//     $("#offListBody").html("");

//     let begin = new Date();
//     let spArr = new Array();
//     let misNum = $("#misNum").val();
//     let spacerInfo = new Array();

//     $.each($("input[name='boxs']:checked"), function () {
//         let index = $(this).val();
//         let tempInfo = new Object();
//         let tmpindex = index - parseInt((index - 1) / 10) * 10;
//         tempInfo.index = index;
//         tempInfo.spacer = $("#resultTable tr:eq(" + tmpindex + ") td:nth-child(6)").text();
//         spArr[index - 1] = tempInfo.spacer;
//         tempInfo.PAM = $("#PAMValue").val();
//         if ($("#spacerLocate").val() === "3' of PAM") {
//             tempInfo.updown = "-";
//         } else {
//             tempInfo.updown = "+";
//         }
//         spacerInfo.push(tempInfo);
//     });

//     let reg = getRegPAMOff($("#PAMValue").val());

//     let subSeq = $("#cdsSequence").val();

//     if (typeof (Worker) !== "undefined") {
//         myWorker = new Worker('js/align.js');
//         $("#myModal").modal("show");
//         myWorker.postMessage([fileTemp, spacerInfo, misNum, reg, subSeq]);
//         myWorker.onmessage = function (e) {
//             if (e.data[0] == 1) {
//                 $("#" + e.data[1] + "Bar").attr("style", "width: " + e.data[2] + "%;");
//             } else if (e.data[0] == 2) {
//                 $("#" + e.data[1] + "Info").text(e.data[2]);
//             } else if (e.data[0] == 3) {
//                 closeModal();
//                 let end = new Date();
//                 let time = end - begin;
//                 console.log("Time consumption is " + time + "ms");
//             } else {
//                 generaOffTableSet(e.data, spArr, "offListBody");
//             }
//         }
//     } else {
//         showSomething("warnNotSupport");
//     }

// }


function generaOffTableSet(InfoSet, spArr, tablename) {
    let $subTable = $("#offResultTable2");
    $subTable.bootstrapTable({
        showColumns: true,
        clickToSelect: true,
        showExport: true,
        exportDataType: "all",
        fixedColumns: true,
        pagination: true,
        sidePagination: "client",
        pageNumber: 1,
        pageSize: 10,
        detailView: true,
        columns: [{
                field: 'index',
                title: 'ID',
                width: 40,
                sortable: true
            }, {
                field: 'spacer',
                title: 'Spacer',
                width: 270,
                sortable: true
            }, {
                field: 'offnum',
                title: 'Off-target',
                width: 112,
                sortable: true
            }

        ],
        onExpandRow: function (index, row, $detail) {
            if (row.subInfo) {
                genSubOffTable($detail, row);
            }
        }
    });
    // debugger;
    let subcontain = new Array();
    let infoArr = new Array();
    for (let tmp2 of InfoSet) {
        let info = tmp2.split("BOSSA");
        if (typeof (infoArr[info[1]]) == "undefined") {
            infoArr[info[1]] = "";
        }
        infoArr[info[1]] += tmp2 + "ATLASBIOINFO";
    }
    for (let i = 1; i < infoArr.length; i++) {
        let sp = infoArr[i].split("ATLASBIOINFO");
        subcontain.push(generaOffTable(sp, spArr[i - 1]));
    }
    $subTable.bootstrapTable('load', subcontain);
}

function genSubOffTable($detail, row) {
    let $subTable = $detail.html('<table></table>').find('table');
    $subTable.bootstrapTable({
        showColumns: true,
        clickToSelect: true,
        showExport: true,
        exportDataType: "all",
        columns: [{
            field: 'miss',
            title: 'Index',
            width: 100,
            sortable: true
        }, {
            field: 'fastatitle',
            title: 'Fasta',
            sortable: true
        }, {
            field: 'strand',
            title: 'Strand',
            width: 80,
            sortable: true
        }, {
            field: 'pos',
            title: 'Locate',
            sortable: true,
            width: 90,
        }, {
            field: 'seq',
            title: 'Sequence',
            width: 280,
            sortable: true
        }, {
            field: 'mismatch',
            title: 'Mismatch',
            sortable: true
        }]
    });
    $subTable.bootstrapTable('load', row.subInfo);
}

function generaOffTable(sp, spacerSeq) {
    // debugger;
    let subInfo = new Array();
    let offResultmain = sp[0].split("BOSSA");
    let index = parseInt(offResultmain[1]);
    if (offResultmain[0] == 0 && sp.length < 3) {
        return {
            'index': index,
            'spacer': spacerSeq,
            'offnum': "<strong style='color:#00DD00'>None</strong>"
        }
    } else {
        let subLable = 1;
        for (let spp in sp) {
            let offResult = sp[spp].split("BOSSA");
            if (offResult.length < 3) {
                continue
            };
            let updown = offResult[3];
            offResult[4] = offResult[4].substr(0, offResult[4].length - 1);
            let coloredSeq = "";
            if ((offResult[6] == "+" && updown == "+") || (offResult[6] == "-" && updown == "-")) {
                offResult[5] = offResult[5].substr(0, offResult[5].length - 9) + "[" + offResult[5].substr(-9, 9) + "]";
            } else {
                offResult[5] = "[" + offResult[5].substr(0, 9) + "]" + offResult[5].substr(9);
            }
            let ntArr = offResult[5].split("");
            for (let j = 0; j < ntArr.length; j++) {
                if (ntArr[j] == "[" || ntArr[j] == "]") {
                    coloredSeq += ntArr[j];
                    continue;
                }
                if (/[ATCG]/.test(ntArr[j])) {
                    coloredSeq += '<strong><font color="red">' + ntArr[j] + '</strong></font>';
                } else {
                    coloredSeq += ntArr[j];
                }
            }
            if ((offResult[6] == "+" && updown == "+") || (offResult[6] == "-" && updown == "-")) {
                coloredSeq += '<span style="background-color: rgba(72,118,255,0.6);">' + offResult[9].toLowerCase() + '</span>';
            } else {
                coloredSeq = '<span style="background-color: rgba(72,118,255,0.6);">' + offResult[9].toLowerCase() + '</span>' + coloredSeq;
            }
            if (offResult[4] == "") {
                offResult[4] = "Null";
            }
            // debugger;
            let spacerPos = offResult[8] + '-' + (parseInt(offResult[8]) + offResult[5].length - 3);
            subInfo.push({
                'miss': subLable,
                'fastatitle': offResult[7],
                'strand': offResult[6],
                'pos': spacerPos,
                'seq': coloredSeq,
                'mismatch': offResult[4]
            });
            subLable++;
        }
        return {
            'index': index,
            'spacer': spacerSeq,
            'offnum': subInfo.length,
            'subInfo': subInfo
        }
    }
}

function generaOffTableNull(offResult, spacerSeq) {
    let index = parseInt(offResult[1]);
    return {
        'index': index,
        'spacer': spacerSeq,
        'offnum': 'None'
    }
}

function disModal() {
    let conCanPro = confirm("Are you sure?");
    if (conCanPro) {
        myWorker.terminate();
        closeModal();
    }
}

function closeModal() {
    if (typeof (myWorker) == undefined) {
        myWorker.terminate();
    }
    $("#myModal").modal("hide");
}

function judgeCheckBox() {
    let count = 0;
    $.each($("input[name='boxs']:checked"), function () {
        count++;
    });
    if (count == 0) {
        let argu = confirm("You didn't choose any Spacer sequence, predict all？");
        if (!argu) {
            return 0;
        } else {
            document.getElementById("allboxs").checked = true;
            $("input[name='boxs']").each(function () {
                $(this).attr("checked", "true");
            });
            return 1;
        }
    } else {
        return 1;
    }
}

function judgeFileType(file, name) {
    hideSomething("warnWrongZip" + name);
    hideSomething("warnWrongFile" + name);
    hideSomething("warnText" + name);
    hideSomething("warnNotSupport" + name);
    let pattZip = /(\.zip)|(\.gz)|(\.7z)|(\.rar)/i;
    let pattAppl = /application/i;
    let pattText = /text.*plain/i;
    if (pattZip.test(file.name)) {
        showSomething("warnWrongZip" + name);
    } else {
        if (/(\.fa)|(\.fasta)/.test(file.name)) {
            showSomething("warnFasta" + name);
            return;
        }
        if (pattAppl.test(file.type)) {
            showSomething("warnWrongFile" + name);
        } else {
            if (pattText.test(file.type) || file.type === "") {
                showSomething("warnText" + name);
            } else {
                showSomething("warnWrongFile" + name);
            }
        }
    }

}

function inputTable() {
    // inputBaseInfo();
    inputTimeCost();
}

function inputTimeCost() {
    $("#timeTableBody").html("");
    generaTR("timeTableBody", ["E.coli CDS", "4.66MB", "0m1s±0.01s", "0m1s±0.02s", "0m0s±0.01s", "0m0s±0.02s"]);
    generaTR("timeTableBody", ["E.coli DNA", "4.50MB", "0m1s±0.01s", "0m0s±0.02s", "0m0s±0.02s", "0m0s±0.01s"]);
    generaTR("timeTableBody", ["S. cerevisiae CDS", "11.10MB", "0m1s±0.04s", "0m1s±0.03s", "0m1s±0.03s", "0m0s±0.01s"]);
    generaTR("timeTableBody", ["S. cerevisiae DNA", "11.79MB", "0m1s±0.04s", "0m1s±0.03s", "0m1s±0.03s", "0m0s±0.01s"]);
    generaTR("timeTableBody", ["C.elegans CDS", "52.00MB", "0m6s±0.20s", "0m4s±0.13s", "0m3s±0.19s", "0m2s±0.01s"]);
    generaTR("timeTableBody", ["D.rerio CDS", "90.86MB", "0m11s±0.23s", "0m7s±0.04s", "0m5s±0.21s", "0m3s±0.01s"]);
    generaTR("timeTableBody", ["C.elegans DNA", "97.24MB", "0m12s±0.24s", "0m8s±0.24s", "0m6s±0.26s", "0m4s±0.02s"]);
    generaTR("timeTableBody", ["M.musculus CDS", "98.64MB", "0m12s±0.41s", "0m8s±0.25s", "0m6s±0.25s", "0m3s±0.01s"]);
    generaTR("timeTableBody", ["A.thaliana CDS", "115.57MB", "0m13s±0.32s", "0m9s±0.28s", "0m6s±0.25s", "0m4s±0.01s"]);
    generaTR("timeTableBody", ["H.sapiens CDS", "146.49MB", "0m18s±0.01s", "0m12s±0.30s", "0m8s±0.36s", "0m5s±0.03s"]);
    generaTR("timeTableBody", ["D.rerio DNA", "1304.19MB", "2m35s±3.57s", "1m46s±3.12s", "1m13s±0.28s", "0m47s±0.25s"]);
    generaTR("timeTableBody", ["M.musculus DNA", "2642.60MB", "5m7s±7.03s", "3m33s±6.42s", "2m32s±6.24s", "1m35s±0.05s"]);
    generaTR("timeTableBody", ["H.sapiens DNA", "2994.31MB", "5m52s±7.89s", "4m2s±7.02s", "2m52s±7.04s", "1m48s±0.20s"]);

}

function generaTR(tablename, infoArr) {
    let paraTemp = $("<tr></tr>");
    for (let i = 0; i < infoArr.length; i++) {
        paraTemp.append('<td>' + infoArr[i] + '</td>');
    }
    paraTemp.appendTo('#' + tablename);
}

// function getRegPAMOff(PAM) {
//     let temp = new Object();
//     if (PAM === "NGG") {
//         temp.regPAM = new RegExp(/\wGG/ig);
//         temp.regPAMr = new RegExp(/CC\w/ig);
//         temp.regInPamPos = [2, 3];
//         temp.regInPamPosR = [1, 2];
//     } else if (PAM === "NGAN") {
//         temp.regPAM = new RegExp(/\wGA\w/ig);
//         temp.regPAMr = new RegExp(/\wTC\w/ig);
//         temp.regInPamPos = [2, 3];
//         temp.regInPamPosR = [2, 3];
//     } else if (PAM === "NGCG") {
//         temp.regPAM = new RegExp(/\wGCG/ig);
//         temp.regPAMr = new RegExp(/CGC\w/ig);
//         temp.regInPamPos = [2, 4];
//         temp.regInPamPosR = [1, 3];
//     } else if (PAM === "NNGRRT") {
//         temp.regPAM = new RegExp(/\w\wG[A|G]{2}T/ig);
//         temp.regPAMr = new RegExp(/A[T|C]{2}C\w\w/ig);
//         temp.regInPamPos = [3, 6];
//         temp.regInPamPosR = [1, 4];
//     } else if (PAM === "NNNRRT") {
//         temp.regPAM = new RegExp(/\w{3}[A|G]{2}T/ig);
//         temp.regPAMr = new RegExp(/A[T|C]{2}\w{3}/ig);
//         temp.regInPamPos = [4, 6];
//         temp.regInPamPosR = [1, 3];
//     } else if (PAM === "TTTV") {
//         temp.regPAM = new RegExp(/TTT[G|A|C]/ig);
//         temp.regPAMr = new RegExp(/[C|T|G]AAA/ig);
//         temp.regInPamPos = [1, 4];
//         temp.regInPamPosR = [1, 4];
//     } else if (PAM === "NG") {
//         temp.regPAM = new RegExp(/\wG/ig);
//         temp.regPAMr = new RegExp(/C\w/ig);
//         temp.regInPamPos = [2, 2];
//         temp.regInPamPosR = [1, 1];
//     }
//     return temp;
// }