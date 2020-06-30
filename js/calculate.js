/*
Apache License Version 2.0

Copyright 2019 Yu Haopeng atlasbioin4@gmail.com atlas_hp@163.com

Licensed under the Apache License, Version 2.0 (the "License");

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

let xlabpos = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100];
let regionChange = function () {
    "use strict";
    let regArr = regionSlid.bootstrapSlider('getValue');
    $('#searchRegionGene').val(regArr[0] + "%-" + regArr[1] + "%");

};

let spacerChange = function () {
    "use strict";
    let spaArr = spacerSlid.bootstrapSlider('getValue');
    let maxValue = Number($("#spacerLenValue1").val());
    if (spaArr[1] > maxValue) {
        $('#spacerValue1').val(spaArr[0] + "-" + maxValue);
    } else {
        $('#spacerValue1').val(spaArr[0] + "-" + spaArr[1]);
    }
    selectModelsChange();
};

let spacerChange2 = function () {
    "use strict";
    let spaArr = spacerSlid2.bootstrapSlider('getValue');
    let maxValue = Number($("#spacerLenValue2").val());
    if (spaArr[1] > maxValue) {
        $('#spacerValue2').val(spaArr[0] + "-" + maxValue);
    } else {
        $('#spacerValue2').val(spaArr[0] + "-" + spaArr[1]);
    }
    selectModelsChange2();
};

let spacerLenChange = function () {

    let spaLenArr = spacerLenSlid.bootstrapSlider('getValue');
    let spaArr = spacerSlid.bootstrapSlider('getValue');
    if (spaArr[1] > spaLenArr) {
        $('#spacerValue1').val(spaArr[0] + "-" + spaLenArr);
    }
    $('#spacerLenValue1').val(spaLenArr);
    selectModelsChange();
};

let ORFLenChange = function () {
    let ORFLenArr = ORFLenSlid.bootstrapSlider('getValue');
    $('#ORFLenValue').val(ORFLenArr);
};

let spacerLenChange2 = function () {
    let spaLenArr = spacerLenSlid2.bootstrapSlider('getValue');
    let spaArr = spacerSlid2.bootstrapSlider('getValue');
    if (spaArr[1] > spaLenArr) {
        $('#spacerValue2').val(spaArr[0] + "-" + spaLenArr);
    }
    $('#spacerLenValue2').val(spaLenArr);
    selectModelsChange2();
};

function calSpacer() {
    "use strict";

    let orfSeq = $("#orfInfoTable tr:eq(7) td:eq(1)").text();
    $("#preResultDiv").show();
    $("#calResultDev").show();

    $("#resWarn").html("");
    let infoDraw = getValueForDraw();

    let searchRegion = $("#searchRegionGene").val().split("-");
    searchRegion[0] = Number(searchRegion[0].split("%")[0]) / 100;
    searchRegion[1] = Number(searchRegion[1].split("%")[0]) / 100;

    let temp = new RegExp(getRegPAM(infoDraw[1]).tpam, "ig");
    let pamPos = new Array();
    let pamPosR = new Array();

    while ((temp.exec(orfSeq) !== null)) {
        pamPos.push(temp.lastIndex);
        temp.lastIndex = temp.lastIndex - infoDraw[1].length + 1;
    }

    let orfSeqM = transStrand(orfSeq);
    while ((temp.exec(orfSeqM) !== null)) {
        pamPosR.push(temp.lastIndex);
        temp.lastIndex = temp.lastIndex - infoDraw[1].length + 1;
    }

    let endReg = new RegExp(/cag|caa|cga/, 'i');
    let endRegm = new RegExp(/cca/, 'i');
    let plusStrand = calPlusStrand(orfSeq, infoDraw, pamPos, endReg, "Plus");
    let minusStrand = calPlusStrand(orfSeqM, infoDraw, pamPosR, endRegm, "Minus");
    debugger;
    drawTwoPie(pamPos.length, plusStrand.length, pamPosR.length, minusStrand.length);
    staticPAMandDraw(plusStrand, minusStrand);
    let result = plusStrand.concat(minusStrand);
    //Add intron calculate
    let addRes = new Array();
    if ($("input[name='intron']:checked").val() == 'cusintron') {
        let inputIntron = getInputCDSRegion();
        result = addIntronPos(result, inputIntron);
        let taddRes = addIntronRes(inputIntron);
        for (let tt in taddRes.addSeq) {
            let tpamPos = new Array();
            let tpamPosR = new Array();
            while ((temp.exec(taddRes.addSeq[tt]) !== null)) {
                tpamPos.push(temp.lastIndex);
                temp.lastIndex = temp.lastIndex - infoDraw[1].length + 1;
            }
            let torfSeqM = transStrand(taddRes.addSeq[tt]);
            while ((temp.exec(torfSeqM) !== null)) {
                tpamPosR.push(temp.lastIndex);
                temp.lastIndex = temp.lastIndex - infoDraw[1].length + 1;
            }
            let tplusStrand = calPlusStrand(taddRes.addSeq[tt], infoDraw, tpamPos, endReg, "Plus");
            let tminusStrand = calPlusStrand(torfSeqM, infoDraw, tpamPosR, endRegm, "Minus");
            let ttres = tplusStrand.concat(tminusStrand);
            if (ttres.length > 0) {
                addRes = addRes.concat(adjustAddPos(ttres, taddRes.addPos[tt], taddRes.len, inputIntron[inputIntron.length - 1]));
            }
        }
    }
    if (addRes.length > 0) {
        result = result.concat(addRes);
    }
    result.sort(function (x, y) {
        return x[10] - y[10];
    });
    addList(result, searchRegion, infoDraw);
    $("#helpnavdiv").affix();

}

function adjustAddPos(addRes, addPos, clen, fastalen) {
    let newAddRes = new Array();
    for (let ttt in addRes) {
        let spPAM = addRes[ttt][1].split("-").map(Number);
        let spSpacer = addRes[ttt][4].split("-").map(Number);
        let spEW = addRes[ttt][6].split("-").map(Number);
        let jumpLab = 1;
        let spPos = parseInt(addRes[ttt][8]);
        let ttarray = [spPAM[0], spPAM[1], spSpacer[0], spSpacer[1], spEW[0], spEW[1], spPos];
        for (let tttt in ttarray) {
            if (ttarray[tttt] - clen > 0) {
                jumpLab = 0;
                ttarray[tttt] = ttarray[tttt] - clen + addPos;
            } else {
                ttarray[tttt] = addPos - clen + ttarray[tttt];
            }

        }
        if (jumpLab) {
            continue;
        }
        addRes[1] = ttarray[0] + "-" + ttarray[1];
        addRes[4] = ttarray[2] + "-" + ttarray[3];
        addRes[6] = ttarray[4] + "-" + ttarray[5];
        addRes[8] = ttarray[6];
        addRes[10] = ttarray[6] / fastalen;
        newAddRes.push(addRes[ttt]);
    }
    return newAddRes;
}

function addIntronRes(intronArr) {
    let tseq = getCDS();
    let infoDraw = getValueForDraw();
    let howmanyCodon = (Math.round(infoDraw[0] / 3) + 1) * 3;
    let tnewSeq = new Array();
    let tnewPos = new Array();
    for (let i = 1; i < intronArr.length - 1; i++) {
        if (i % 2 == 0) {
            tnewSeq.push(tseq.seq[0].substring(intronArr[i] - 1 - howmanyCodon, intronArr[i] + howmanyCodon - 1));
            tnewPos.push(intronArr[i] - 1);
        } else {
            tnewSeq.push(tseq.seq[0].substring(intronArr[i] - howmanyCodon, intronArr[i] + howmanyCodon));
            tnewPos.push(intronArr[i]);
        }
    }
    return {
        'addSeq': tnewSeq,
        'addPos': tnewPos,
        "len": howmanyCodon
    };
}

function getIntronPos(intronArr) {
    let newPosArr = new Array();
    for (let i = 2; i < intronArr.length; i += 2) {
        newPosArr.push(intronArr[i - 1]);
        //这里我考虑过，1..2,5..6，5-2-1=2.
        newPosArr.push(intronArr[i] - intronArr[i - 1] - 1);
    }
    return newPosArr;
}

function addIntronPos(result, cdsArr) {
    let posArr = getIntronPos(cdsArr);
    let newRes = new Array();
    for (let i in result) {
        let spPAM = result[i][1].split("-").map(Number);
        let spSpacer = result[i][4].split("-").map(Number);
        let jumpLab = 0;
        for (let k in cdsArr) {
            if ((spPAM[0] <= cdsArr[k] && spPAM[1] >= cdsArr[k])) {
                jumpLab = 1;
                break;
            }
            if ((spSpacer[0] <= cdsArr[k] && spSpacer[1] >= cdsArr[k])) {
                jumpLab = 1;
                break;
            }
        }
        if (jumpLab) {
            continue;
        }
        let spEW = result[i][6].split("-").map(Number);
        let teditpos = parseInt(result[i][8]);
        for (let j = 0; j < posArr.length; j += 2) {
            if (spPAM[0] > posArr[j]) {
                spPAM[0] += posArr[j + 1];
            }
            if (spPAM[1] > posArr[j]) {
                spPAM[1] += posArr[j + 1];
            }
            if (spEW[0] > posArr[j]) {
                spEW[0] += posArr[j + 1];
            }
            if (spEW[1] > posArr[j]) {
                spEW[1] += posArr[j + 1];
            }
            if (spSpacer[0] > posArr[j]) {
                spSpacer[0] += posArr[j + 1];
            }
            if (spSpacer[1] > posArr[j]) {
                spSpacer[1] += posArr[j + 1];
            }
            if (teditpos > posArr[j]) {
                teditpos += posArr[j + 1];
            }
        }
        result[i][1] = spPAM.join("-");
        result[i][4] = spSpacer.join("-");
        result[i][6] = spEW.join("-");
        result[i][8] = teditpos;
        result[i][10] = teditpos / (cdsArr[cdsArr.length - 1] - cdsArr[0] + 1);
        newRes.push(result[i]);
    }
    return newRes;
}

function getInputCDSRegion() {
    let intronArr = new Array();
    let tintrontext = $("#cusintrontext").val();
    tintrontext = tintrontext.replace(/\s|\r|\t|\n/, "");
    let rawArr = tintrontext.split(",");
    for (let i in rawArr) {
        let tarrarr = rawArr[i].split("..")
        intronArr.push(parseInt(tarrarr[0]));
        intronArr.push(parseInt(tarrarr[1]));
    }

    return intronArr;
}

function staticPAMandDraw(plusStrand, minusStrand) {
    let tname = new Array();
    let tpos = new Array();
    tname.push('Plus strand');
    tname.push('Minus strand');
    let ptpos = Array.apply(null, Array(100)).map(function (v, i) {
        return 0;
    });
    let mtpos = Array.apply(null, Array(100)).map(function (v, i) {
        return 0;
    });
    for (let ptmp in plusStrand) {
        let pos = plusStrand[ptmp][10];
        pos = parseInt(pos * 100);
        ptpos[pos]++;
    }
    for (let mtmp in minusStrand) {
        let pos = minusStrand[mtmp][10];
        pos = parseInt(pos * 100);
        mtpos[pos]++;
    }
    tpos.push(ptpos);
    tpos.push(mtpos);
    let data = {
        'name': tname,
        'count': tpos
    };
    drawPAMStaBar(data);
}

function getCDS() {
    let seq = $("#cdsSequence").val();
    let seqArr = seq.split("\n");
    let cds = new Array();
    let cdsName = new Array();
    let j = 0;
    if (seqArr[0].match(/^>/)) {
        cdsName[j] = seqArr[0];
    } else {
        cds[0] = seqArr[0].toUpperCase();
    }
    for (let i = 1; i < seqArr.length; i++) {
        if (seqArr[i].match(/^>/)) {
            j++;
            cdsName[j] = seqArr[i];
            continue;
        }
        seqArr[i] = seqArr[i].replace(/\s|\r|\t|\n/, "").toUpperCase();
        if (cds[j]) {
            cds[j] = cds[j] + seqArr[i];
        } else {
            cds[j] = seqArr[i];
        }
    }
    return {
        "seq": cds,
        "seqName": cdsName
    };
}

function getRegPAM(PAM) {
    let tpam = PAM.split("");
    for (let t in tpam) {
        if (/[A,T,C,G]/i.test(tpam[t])) {
            continue;
        }
        if (/N/i.test(tpam[t])) {
            tpam[t] = '\\w';
        } else if ((/R/i.test(tpam[t]))) {
            tpam[t] = "[AGR]";
        } else if ((/Y/i.test(tpam[t]))) {
            tpam[t] = "[CTY]";
        } else if ((/M/i.test(tpam[t]))) {
            tpam[t] = "[ACM]";
        } else if ((/K/i.test(tpam[t]))) {
            tpam[t] = "[GTK]";
        } else if ((/S/i.test(tpam[t]))) {
            tpam[t] = "[GCS]";
        } else if ((/W/i.test(tpam[t]))) {
            tpam[t] = "[ATW]";
        } else if ((/H/i.test(tpam[t]))) {
            tpam[t] = "[ATCWMYH]";
        } else if ((/B/i.test(tpam[t]))) {
            tpam[t] = "[GTCKSYB]";
        } else if ((/V/i.test(tpam[t]))) {
            tpam[t] = "[GACRSMV]";
        } else if ((/D/i.test(tpam[t]))) {
            tpam[t] = "[GATRKWD]";
        }
    }
    let tpamr = transStrand(PAM).split("");
    for (let tr in tpamr) {
        if (/[A,T,C,G]/i.test(tpamr[tr])) {
            continue;
        }
        if (/N/i.test(tpamr[tr])) {
            tpamr[tr] = '\\w';
        } else if ((/R/i.test(tpamr[tr]))) {
            tpamr[tr] = "[AGR]";
        } else if ((/Y/i.test(tpamr[tr]))) {
            tpamr[tr] = "[CTY]";
        } else if ((/M/i.test(tpamr[tr]))) {
            tpamr[tr] = "[ACM]";
        } else if ((/K/i.test(tpamr[tr]))) {
            tpamr[tr] = "[GTK]";
        } else if ((/S/i.test(tpamr[tr]))) {
            tpamr[tr] = "[GCS]";
        } else if ((/W/i.test(tpamr[tr]))) {
            tpamr[tr] = "[ATW]";
        } else if ((/H/i.test(tpamr[tr]))) {
            tpamr[tr] = "[ATCWMYH]";
        } else if ((/B/i.test(tpamr[tr]))) {
            tpamr[tr] = "[GTCKSYB]";
        } else if ((/V/i.test(tpamr[tr]))) {
            tpamr[tr] = "[GACRSMV]";
        } else if ((/D/i.test(tpamr[tr]))) {
            tpamr[tr] = "[GATRKWD]";
        }
    }
    // let treg = new RegExp(tpam.join(""), "ig");
    return {
        'tpam': tpam.join(""),
        'tpamr': tpamr.join(""),
        'len': PAM.length
    };
}

function showHelp() {
    $('#cbeiTab a[href="#Help"]').tab('show');
    inputTable();
    scrollTo(0, 0);
}


function resultDis(bestResult, len, screenRes, searchRegion, infoDraw) {
    "use strict";
    $("#beDiscript").html("");
    $("#beDiscript").append($('<li><strong style="color:red">' + len + '</strong> potential editing sites, <strong style="color:red">' + screenRes + '</strong> in the search area (' + Math.round(searchRegion[0] * 10000) / 100 + '%' + '-' + Math.round(searchRegion[1] * 10000) / 100 + '%' + ').</li>'));
    $("#sequenceDiscript").show();
    $("#calResultDev").show();
    if (typeof (bestResult) == 'undefined' && len == 0) {
        $("#beDiscript").append($('<li><strong style="color:red" >Oh, CRISPR-CBEI has not found a potential site, please check the input sequence.</strong></li>'));
        $("#sequenceDiscript").hide();
        $("#calResultDev").hide();
        return 0;
    }
    if (len > 0 && screenRes == 0) {
        $("#beDiscript").append($('<li><strong style="color:blue">CRISPR-CBEI found potential sites, but outside the search area you set up. Please try to reset the search area.</strong></li>'));
    } else {
        $("#beDiscript").append($('<li><strong>The closest potential editing site to the 5\' end is:</strong></li>'));
        drawPAM("resPAM", infoDraw[0], infoDraw[1], infoDraw[2], infoDraw[3], [infoDraw[4], infoDraw[5]], infoDraw[6]);
        let $bestTableLi = $('<li></li>');
        let $bestTable = $bestTableLi.html('<table></table>').find('table');
        $bestTable.bootstrapTable({
            showColumns: true,
            clickToSelect: true,
            showExport: true,
            exportDataType: "all",
            columns: [{
                field: 'item',
                title: 'Items',
                sortable: true
            }, {
                field: 'value',
                title: 'Value',
                sortable: true
            }]
        });
        let subcontain = new Array();
        subcontain.push({
            'item': 'Strand',
            'value': bestResult[2]
        });
        subcontain.push({
            'item': 'Spacer details',
            'value': bestResult[3]
        });
        subcontain.push({
            'item': 'Spacer region',
            'value': bestResult[4]
        });
        subcontain.push({
            'item': 'Edit window',
            'value': bestResult[5]
        });
        subcontain.push({
            'item': 'Edit region',
            'value': bestResult[6]
        });
        // subcontain.push({
        //     'item': 'Edit codon',
        //     'value': bestResult[7]
        // });
        subcontain.push({
            'item': 'Edit position',
            'value': bestResult[8]
        });
        subcontain.push({
            'item': 'Patten',
            'value': bestResult[9]
        });
        subcontain.push({
            'item': 'Location',
            'value': bestResult[10]
        });
        $bestTable.bootstrapTable('load', subcontain);
        $("#beDiscript").append($bestTableLi);
        return 1;
    }
}

function resetAll() {
    $("#ORFidentificationDiv").hide();
    $("#CBEIpredictionDIV").hide();
    $("#preResultDiv").hide();
    $("#calResultDev").hide();
    $("#childAdORF").hide();
    $("#childAd").hide();
    $("#cdsSequence").val("");
    $("#spacerValue").val("4-8");
    $("#searchRegionGene").val("0%-50%");
    $("#ORFLenValue").val("75");
    $("input[value='ORFoption1']").prop('checked', true);
    $("input[value='AUG']").prop('checked', true);
    $("input[value='sel']").prop('checked', true);
    $("#ORFminlength").bootstrapSlider("setValue", 75);
    $("#spacerSlider").bootstrapSlider("setValue", [4, 8]);
    $("#searchRegionSlider").bootstrapSlider("setValue", [0, 50]);
    $("#helpnavdiv").affix();
}

function resetCBEI() {
    $("#preResultDiv").hide();
    $("#calResultDev").hide();
    $("#spacerValue").val("4-8");
    $("input[value='sel']").prop('checked', true);
    $("#PAMmodels1").selectpicker('val', 'BE:"NGG,Spacer:20nt,Edit:[4-8]"');
    $("#childAd").hide();
    $("#searchRegionGene").val("0%-50%");
    $("#ORFminlength").bootstrapSlider("setValue", 75);
    $("#spacerSlider").bootstrapSlider("setValue", [4, 8]);
    $("#searchRegionSlider").bootstrapSlider("setValue", [0, 50]);
    selectModel(1);
}

function resetInput() {
    $("#cdsSequence").val("");
    $("#ORFLenValue").val("75");
    $("input[value='ORFoption1']").prop('checked', true);
    $("input[value='AUG']").prop('checked', true);
    $("#childAdORF").hide();

}

function returnToSub() {
    "use strict";
    $('#cbeiTab a[href="#submitCDS"]').tab('show');
    scrollTo(0, 0);
    // hideSomething("warnNoResult");
    hideSomething("desirPAMDiscript");
    hideSomething("stopCodonDiscript");
    // hideSomething("warnGroup");
    hideSomething("calResultDev");
    // hideSomething("warnSuccess");
    hideSomething("stopCodonDiscript");
    // hideSomething("warnNoLastStop");
    hideSomething("stopCodonDev");
    // hideSomething("warnMutiStop");
}

function showSomething(name) {
    "use strict";
    document.getElementById(name).style.display = "inline";
    document.getElementById(name).style.visibility = "visible";
}

function hideSomething(name) {
    "use strict";
    document.getElementById(name).style.display = "none";
    document.getElementById(name).style.visibility = "hidden";
}

function isNotInSearchRegion(value, searchRegion) {
    "use strict";
    if (value > searchRegion[0] && value < searchRegion[1]) {
        return false;
    } else {
        return true;
    }
}

function resultOut(id, res) {
    return {
        "id": id,
        "strand": res[2],
        'pam': res[0],
        'pamr': res[1],
        'spacer': res[11],
        'spacerpos': res[4],
        'region': res[5],
        'regionpos': res[6],
        'codon': res[7],
        'editp': res[8],
        'pattern': res[9],
        'location': res[10],
        'nolabspacer': res[3],
        'purespacer': res[12]
    };
}

function addList(result, searchRegion, infoDraw) {
    "use strict";
    let orilength = result.length;
    result = trimmedResultFromRegion(result, searchRegion);
    let tResult = new Array();
    for (let re in result) {
        let tre = parseInt(re) + 1;
        tResult.push(resultOut(tre, result[re]));
    }
    $("#resultTable").bootstrapTable("load", tResult);

    let signal = resultDis(result[0], orilength, result.length, searchRegion, infoDraw);
    if (signal == 0) {

    }
    return result.length;
}

function trimmedResultFromRegion(result, searchRegion) {
    let newArr = new Array();
    for (let sp = 0; sp < result.length; sp++) {
        if (isNotInSearchRegion(result[sp][10], searchRegion)) {
            continue;
        } else {
            result[sp][10] = Math.round(result[sp][10] * 10000) / 100 + "%";
            newArr.push(result[sp]);
        }
    }
    return newArr;
}

function genResultTableHeader() {
    let $rtable = $("#resultTable");
    $rtable.bootstrapTable({
        striped: true,
        pagination: true,
        sidePagination: "client",
        toolbar: "#resultToolbar",
        pageNumber: 1,
        pageSize: 10,
        // strictSearch: true,
        showColumns: true,
        clickToSelect: true,
        uniqueId: "id",
        showExport: true,
        exportDataType: "all",
        detailView: true,
        columns: [{
                field: 'cbox',
                checkbox: 'true'
            },
            {
                field: 'id',
                title: 'ID',
                sortable: true,
                align: "center",
                width: "50"
            }, {
                field: 'strand',
                title: 'Strand',
                sortable: true,
                align: "center",
                width: "88"
            }, {
                field: 'spacer',
                title: 'Spacer',
                sortable: true,
                align: "center",
                width: "420"
            }, {
                field: 'pattern',
                title: 'Pattern',
                align: "center",
                sortable: true
            }, {
                field: 'location',
                title: 'Location',
                align: "center",
                sortable: true
            }
        ],
        onExpandRow: function (index, row, $detail) {
            genSubResTable($detail, row);

        },
    });

    $(function () {
        $("#ebut1").click(function () {
            $rtable.bootstrapTable('expandAllRows');
        });
        $("#cbut1").click(function () {
            $rtable.bootstrapTable('collapseAllRows')
        });
    });
}

function genSubResTable($detail, row) {
    let $subTable = $detail.html('<table></table>').find('table');
    $subTable.bootstrapTable({
        showColumns: true,
        clickToSelect: true,
        showExport: true,
        exportDataType: "all",
        columns: [{
            field: 'item',
            title: 'Items',
            sortable: true
        }, {
            field: 'value',
            title: 'Value',
            sortable: true
        }]
    });
    let subcontain = new Array();
    subcontain.push({
        'item': 'Strand',
        'value': row.strand
    });
    subcontain.push({
        'item': 'Spacer details',
        'value': row.nolabspacer
    });
    subcontain.push({
        'item': 'Spacer region',
        'value': row.spacerpos
    });
    subcontain.push({
        'item': 'Edit window',
        'value': row.region
    });
    subcontain.push({
        'item': 'Edit region',
        'value': row.regionpos
    });
    // subcontain.push({
    //     'item': 'Edited codon',
    //     'value': row.codon
    // });
    subcontain.push({
        'item': 'Edited position',
        'value': row.editp
    });
    subcontain.push({
        'item': 'Edit patten',
        'value': row.pattern
    });
    subcontain.push({
        'item': 'PAM',
        'value': row.pam
    });
    subcontain.push({
        'item': 'PAM region',
        'value': row.pamr
    });
    subcontain.push({
        'item': 'Location',
        'value': row.location
    });
    $subTable.bootstrapTable('load', subcontain);
}

function allcheck() {
    let nn = $("#allboxs").is(":checked");
    if (nn == true) {
        let namebox = $("input[name^='boxs']");
        for (let i = 0; i < namebox.length; i++) {
            namebox[i].checked = true;
        }
    }
    if (nn == false) {
        let namebox = $("input[name^='boxs']");
        for (let i = 0; i < namebox.length; i++) {
            namebox[i].checked = false;
        }
    }

}

function drawPAM(canV, lenToDraw, pamseqP, pamseqM, pamPos, editRegArr, strand) {
    "use strict";
    // figPAM 23 NGG NCC 21 Array [ 4, 8 ] +
    // figPAM 20 NGG 21 [4,8] +
    let cc = document.getElementById(canV);
    if (cc.getContext) {
        let cxt = cc.getContext("2d");
        cc.height = cc.height;
        cxt.font = "15px Arial";
        biankuang(cxt, lenToDraw);
        draw53(cxt, 30, "5'", "3'");
        cxt.beginPath();
        cxt.lineWidth = 4;
        cxt.moveTo(45, 30);
        cxt.lineTo(45 + (lenToDraw + 1) * 14, 30);
        cxt.stroke();
        cxt.moveTo(45, 55);
        cxt.lineTo(45 + (lenToDraw + 1) * 14, 55);
        cxt.stroke();

        cxt.beginPath();
        cxt.lineWidth = 2;
        for (let i = 1; i <= lenToDraw; i++) {
            cxt.moveTo(45 + i * 14, 30);
            cxt.lineTo(45 + i * 14, 55);
            cxt.stroke();
        }
        draw53(cxt, 45 + (lenToDraw + 1) * 14 + 5, "3'", "5'");
        cxt.globalAlpha = 0.4;
        cxt.fillStyle = "#4876FF";
        cxt.clearRect(45 + pamPos * 14 - 7, 20, (pamseqP.length) * 14, 20);
        cxt.fillRect(45 + pamPos * 14 - 7, 20, (pamseqP.length) * 14, 20);
        cxt.clearRect(45 + pamPos * 14 - 7, 45, (pamseqP.length) * 14, 20);
        cxt.fillRect(45 + pamPos * 14 - 7, 45, (pamseqP.length) * 14, 20);
        cxt.globalAlpha = 1;

        addSpacerNum(cxt, lenToDraw - pamseqP.length, pamseqP.length, strand);
        addEditRegion(cxt, editRegArr, pamseqP.length, strand);

        let tmp = pamseqP.split("");
        pamseqP = tmp.join(" ");
        tmp = pamseqM.split("");
        pamseqM = tmp.join(" ");
        cxt.fillText(pamseqP, 45 + pamPos * 14 - 5, 35);
        cxt.fillText(pamseqM, 45 + pamPos * 14 - 5, 60);

    }
}

function addEditRegion(cxt, editRegArr, pamLen, strand) {
    "use strict";
    if (strand === "+") {
        if (editRegArr[0] === editRegArr[1]) {
            if (editRegArr[0] > 9) {
                cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + 7, 24);
            } else {
                cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + 9, 24);
            }
        } else {
            cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + 7, 24);
            cxt.fillText(editRegArr[1], 45 + (editRegArr[1] - 1) * 14 + 7, 24);
        }
        cxt.globalAlpha = 0.3;
        cxt.fillStyle = "#FF3030";
        if (editRegArr[0] === editRegArr[1]) {
            if (editRegArr[0] > 9) {
                cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 7, 10, 18, 30);
            } else {
                cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 7, 10, (editRegArr[1] - editRegArr[0] + 1) * 14, 30);
            }

        } else {
            cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 7, 10, (editRegArr[1] - editRegArr[0] + 1) * 14, 30);
        }
    } else {
        if (editRegArr[0] === editRegArr[1]) {
            if (editRegArr[0] > 9) {
                cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + pamLen * 14 + 7, 24);
            } else {
                cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + pamLen * 14 + 9, 24);
            }
        } else {
            cxt.fillText(editRegArr[0], 45 + (editRegArr[0] - 1) * 14 + pamLen * 14 + 7, 24);
            cxt.fillText(editRegArr[1], 45 + (editRegArr[1] - 1) * 14 + pamLen * 14 + 7, 24);
        }
        cxt.globalAlpha = 0.3;
        cxt.fillStyle = "#FF3030";
        if (editRegArr[0] === editRegArr[1]) {
            if (editRegArr[0] > 9) {
                cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 2 + pamLen * 14 + 7, 10, 15, 30);
            } else {
                cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 2 + pamLen * 14 + 7, 10, (editRegArr[1] - editRegArr[0] + 1) * 14, 30);
            }

        } else {
            cxt.fillRect(45 + (editRegArr[0] - 1) * 14 + 2 + pamLen * 14 + 7, 10, (editRegArr[1] - editRegArr[0] + 1) * 14, 30);
        }
    }

    cxt.globalAlpha = 1;
    cxt.fillStyle = "#000";
}

function addSpacerNum(cxt, spacerLen, pamLen, strand) {
    "use strict";
    cxt.fillStyle = "#000";
    if (strand === "+") {
        cxt.fillText(1, 45 + 9, 70);
        cxt.fillText(spacerLen, 45 + spacerLen * 14 - 9, 70);
        cxt.font = "bold 15px Arial";
        cxt.fillText("Spacer", 40 + spacerLen * 7, 75);
        cxt.font = "15px Arial";
        cxt.globalAlpha = 0.2;
        cxt.fillStyle = "#66CD00";
        cxt.fillRect(45 + 9, 20, spacerLen * 14 - 4, 60);
    } else {
        cxt.fillText(1, 45 + 9 + pamLen * 14, 70);
        cxt.fillText(spacerLen, 45 + spacerLen * 14 - 9 + pamLen * 14, 70);
        cxt.font = "bold 15px Arial";
        cxt.fillText("Spacer", 40 + spacerLen * 7 + pamLen * 7, 75);
        cxt.font = "15px Arial";
        cxt.globalAlpha = 0.2;
        cxt.fillStyle = "#66CD00";
        cxt.fillRect(45 + 9 + pamLen * 14, 20, spacerLen * 14 - 4, 60);
    }
    cxt.globalAlpha = 1;
    cxt.fillStyle = "#000";
}

function draw53(cxt, colP, name1, name2) {
    "use strict";
    cxt.fillText(name1, colP, 35);
    cxt.fillText(name2, colP, 60);
}

function biankuang(cxt, lenToDraw) {
    "use strict";
    cxt.beginPath();
    cxt.moveTo(20, 0);
    cxt.lineTo(20, 90);
    cxt.lineTo(45 + lenToDraw * 14 + 38, 90);
    cxt.lineTo(45 + lenToDraw * 14 + 38, 0);
    cxt.lineTo(20, 0);
    cxt.stroke();
    cxt.closePath();
}

function transStrand(seq) {
    "use strict";
    let baseT = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "U": "A",
        ",": ",",
        "-": "-",
        "[": "]",
        "]": "[",
        "(": ")",
        ")": "(",
        "N": "N",
        ">": "<",
        "<": ">",
        "R": "Y",
        "Y": "R",
        "M": "K",
        "S": "S",
        "W": "W",
        "H": "D",
        "D": "H",
        "V": "B",
        "B": "V"
    };
    let rSeq = seq.split("").reverse();

    for (let i = 0; i < rSeq.length; i++) {
        if (/[a-z]/.test(rSeq[i])) {
            rSeq[i] = baseT[rSeq[i].toUpperCase()].toLowerCase();
        } else {
            rSeq[i] = baseT[rSeq[i]];
        }
    }
    return rSeq.join("");
}



function calPlusStrand(cds, infoDraw, pamPosLast, endReg, mode) {
    let preStop = new Array();
    let direc = infoDraw[6];
    let spacerLen = infoDraw[0] - infoDraw[1].length;
    for (let k = 0; k < pamPosLast.length; k++) {
        let pamInfo = new Array([11]);
        let b;
        let e;
        let spb;
        let spe;
        if (direc === "+") {
            if (pamPosLast[k] - infoDraw[0] > -1) {
                spb = pamPosLast[k] - infoDraw[0] + 1;
                spe = pamPosLast[k] - infoDraw[1].length;
            } else {
                continue;
            }
            b = pamPosLast[k] - infoDraw[0] + infoDraw[4];
            e = pamPosLast[k] - infoDraw[0] + infoDraw[5];
            pamInfo[1] = cds.substr(pamPosLast[k] - infoDraw[1].length, infoDraw[1].length);
        } else {
            if (pamPosLast[k] + spacerLen > cds.length) {
                continue;
            }
            spb = pamPosLast[k] + 1;
            spe = pamPosLast[k] + spacerLen;
            b = pamPosLast[k] + infoDraw[4];
            e = pamPosLast[k] + infoDraw[5];
            pamInfo[1] = cds.substr(pamPosLast[k] - infoDraw[1].length, infoDraw[1].length);
        }
        if (mode == "Minus") {
            pamInfo[5] = (cds.length - spb + 1) + "-" + (cds.length - spe + 1);
            pamInfo[7] = (cds.length - b + 1) + "-" + (cds.length - e + 1);
            pamInfo[2] = (cds.length - pamPosLast[k] + infoDraw[1].length) + "-" + (cds.length - pamPosLast[k] + 1);
        } else {
            pamInfo[5] = spb + "-" + spe;
            pamInfo[7] = b + "-" + e;
            pamInfo[2] = (pamPosLast[k] - infoDraw[1].length + 1) + "-" + pamPosLast[k];
        }


        pamInfo[3] = infoDraw[6];
        let pureSpacer = cds.substr(spb - 1, spe - spb + 1);

        let orispb = spb;
        let orispe = spe;
        let longspb = spb;
        let longspe = spe;

        if (orispb % 3 === 0) {
            longspb = longspb - 2;
        } else if (orispb % 3 === 2) {
            longspb = longspb - 1;
        }

        if (orispe % 3 == 1) {
            longspe = longspe + 2;
        } else if (orispe % 3 == 2) {
            longspe = longspe + 1;
        }
        let longspacer = cds.substr(longspb - 1, longspe - longspb + 1);
        let longcodon = longspacer.length / 3;
        let tedit = new Array();
        let tcodonall = new Array();
        for (let tsp = 0; tsp < longcodon; tsp++) {
            let tcodon = longspacer.substr(tsp * 3, 3);
            if (tcodon.match(endReg)) {
                tedit.push(longspb + tsp * 3);
            }
            tcodonall.push(tcodon);
        }
        let tcodonallori = tcodonall.concat();
        let endPam;

        if (direc == "+") {
            tcodonall[parseInt((b - longspb) / 3)] = spCodon(tcodonallori[parseInt((b - longspb) / 3)], toThree(b % 3), "b", '<strong>[</strong><span style="background-color: wheat;">');
            tcodonall[parseInt((e - longspb) / 3)] = spCodon(tcodonallori[parseInt((e - longspb) / 3)], toThree(e % 3), "e", "</span><strong>]</strong>");

            tcodonall[0] = spCodon(tcodonall[0], orispb - longspb + 1, "b", '<strong>{</strong><span style="background-color: rgba(102,205,0,0.3);">');
            tcodonall[tcodonall.length - 1] = spCodon(tcodonall[tcodonall.length - 1], 3 - longspe + orispe, "e", '</span><strong>}|</strong><span style="background-color: rgba(72,118,255,0.6);">');
        } else {
            tcodonall[parseInt((b - longspb) / 3)] = spCodon(tcodonallori[parseInt((b - longspb) / 3)], toThree(b % 3), "b", '<strong>[</strong><span style="background-color:wheat;">');
            tcodonall[parseInt((e - longspb) / 3)] = spCodon(tcodonall[parseInt((e - longspb) / 3)], toThree(e % 3), "e", "</span><strong>]</strong>");
            tcodonall[0] = spCodon(tcodonall[0], orispb - longspb + 1, "b", '</span><strong>|{</strong><span style="background-color: rgba(102,205,0,0.3);">');
            tcodonall[tcodonall.length - 1] = spCodon(tcodonall[tcodonall.length - 1], 3 - longspe + orispe, "e", '</span><strong>}</strong>');
        }
        let editPatten = "";
        let tlocation = -1.1;
        let tios = new Array();
        let noskipReg = new RegExp('[ATCG]');
        // debugger;
        for (let ttedit in tedit) {
            let teditCodonIndex = parseInt((tedit[ttedit] - longspb) / 3);

            for (let pedit = 0; pedit < 3; pedit++) {
                let teditp = pedit + tedit[ttedit];
                let spCodon = tcodonallori[teditCodonIndex].split("");
                if (teditp >= b && teditp <= e) {
                    if (/c/i.test(spCodon[pedit])) {
                        tios.push(teditp);
                        tcodonall[teditCodonIndex] = changeLabel(tcodonall[teditCodonIndex], pedit, '<strong style="color: Crimson">(c->t)</strong>', noskipReg)
                        // tcodonall[teditCodonIndex] = changeCondonStrSingle(tcodonall[teditCodonIndex], pedit+1, '<strong style="color: wheat">(C->t)</strong>');
                    }
                }
            }
        }
        if (tios.length == 0) {
            continue;
        }

        let teditpatten = new Array();
        for (let ii in tios) {
            teditpatten.push(cds.substr(tios[ii] - 2, 2));
        }
        editPatten = teditpatten.join(",");

        if (mode == 'Minus') {
            tlocation = (cds.length - tios[0] + 1) / cds.length;
        } else {
            tlocation = tios[0] / cds.length;
        }

        for (let iii in tios) {
            if (mode == 'Minus') {
                tios[iii] = cds.length - tios[iii] + 1;
            }
        }
        let editPos = tios.join(",");

        if (editPos.length == 0) {
            continue;
        }
        let detailSeq = tcodonall.join(",");
        if (direc == "+") {
            endPam = cds.substr(longspe, infoDraw[1].length - longspe + orispe) + "</span>";
            detailSeq = detailSeq + ',' + endPam;
        } else {
            endPam = '<span style="background-color: rgba(72,118,255,0.6);">' + cds.substr(longspb - (infoDraw[1].length + 1 - orispb + longspb), infoDraw[1].length - orispb + longspb);
            detailSeq = endPam + "," + detailSeq;
        }
        // debugger;
        detailSeq = detailSeq.toUpperCase();
        tcodonall[0] = tcodonall[0].replace(/^\w+<strong>/g, "<strong>");
        tcodonall[tcodonall.length - 1] = tcodonall[tcodonall.length - 1].replace(/<\/strong>\w+$/g, "</strong>");
        let nolablespacer = detailSeq.split(",").join("");
        if (/\{(.*?)\}/.test(nolablespacer)) {
            nolablespacer = /\{(.*?)\}/.exec(nolablespacer)[1].trim();
        }


        nolablespacer = nolablespacer.replace(/\{|\}|\[|\]|\|/g, "");
        nolablespacer = nolablespacer.replace(/\(C->T\)/ig, "C");
        nolablespacer = nolablespacer.replace(/10220500\.3/, "102,205,0,0.3");
        let teditseq = "";
        if (/\[(.*?)\]/.test(detailSeq)) {
            teditseq = /\[(.*?)\]/.exec(detailSeq)[1].trim();
        }
        teditseq = teditseq.replace(/\(C->T\)/ig, "C");

        // "strand": pamInfo[3]
        // 'pam': pamInfo[1]
        //  'pamr': pamInfo[2],
        // 'spacer': pamInfo[4]
        // 'spacerpos': pamInfo[5]
        // 'region': pamInfo[6], 
        // 'regionpos': pamInfo[7] , 
        // 'codon': editCodon.join(";")
        // 'editp': editpos.join(";"), 
        // 'pattern'editPatten[tindex]: 
        // 'location': res[10],
        // 'nolabspacer':res[11] };

        let tempArr = new Array(pamInfo[1], pamInfo[2], mode, detailSeq, pamInfo[5], teditseq, pamInfo[7], 'atlas', editPos, editPatten, tlocation, nolablespacer, pureSpacer);
        preStop.push(tempArr);

    }
    return preStop;
}

function changeLabel(seq, addPos, addStr, reg) {
    let count = -1;
    let newSeq = '';
    let errLab = 1;
    let spseq = seq.split("");
    for (let i in spseq) {
        if (reg.test(spseq[i])) {
            count++;
            if (count == addPos) {
                newSeq = newSeq + addStr;
                errLab = 0;
            } else {
                newSeq = newSeq + spseq[i];
            }
        } else {
            newSeq = newSeq + spseq[i];
        }
    }
    if (errLab) {
        newSeq = newSeq + addStr;
    }
    return newSeq;
}

function addLabel(seq, addPos, addStr, reg) {
    let count = -1;
    let newSeq = '';
    let errLab = 1;
    let spseq = seq.split("");
    for (let i in spseq) {
        if (reg.test(spseq[i])) {
            count++;
            if (count == addPos) {
                newSeq = newSeq + addStr + spseq[i];
                errLab = 0;
            } else {
                newSeq = newSeq + spseq[i];
            }
        } else {
            newSeq = newSeq + spseq[i];
        }
    }
    if (errLab) {
        newSeq = newSeq + addStr;
    }
    return newSeq;
}

function toThree(num) {
    if (num == 0) {
        num = 3;
    }
    return num;
}

function spCodon(beEdit, numOfEdit, b, change) {
    let arr = beEdit.split("");
    let index = 0;
    for (let a in arr) {
        if (arr[a].match(/[^A-Z]/)) {
            continue;
        }
        index++;
        if (numOfEdit == index) {
            if (b == "b") {
                arr[a] = change + arr[a];
                break;
            } else if (b == "e") {
                arr[a] = arr[a] + change;
                break;
            }
        }
    }

    return arr.join("");
}

function changeCondonStrSingle(codon, pos, change) {
    let ccodon = codon;
    if (pos == 1) {
        ccodon = change + codon.substr(pos);
    } else if (pos == codon.length) {
        ccodon = codon.substr(0, pos - 1) + change;
    } else {
        ccodon = codon.substr(0, pos - 1) + change + codon.substr(pos);
    }
    return ccodon;
}

function judgeSpacerEnd(tspacer, beg, e) {
    let espa = "None";
    let cspa = "None";
    if (e % 3 === 1) {
        espa = tspacer.substr(tspacer.length - 1, 1);
        e = e - 1;
        cspa = tspacer.substr(beg, tspacer.length - 1 - beg);
    } else if (e % 3 === 2) {
        espa = tspacer.substr(tspacer.length - 2, 2);
        e = e - 2;
        cspa = tspacer.substr(beg, tspacer.length - 2 - beg);
    } else {
        cspa = tspacer.substr(beg, tspacer.length - beg);
    }
    return {
        "espa": espa,
        "cspa": cspa,
        "spe": e
    };
}

function exportListCSV(tablename, fileName) {
    "use strict";
    $('#' + tablename).tableExport({
        type: 'csv',
        fileName: fileName
    });
}

function exportListTSV(tablename, fileName) {
    "use strict";
    $('#' + tablename).tableExport({
        type: 'tsv',
        fileName: fileName
    });
}

function exportListJSON(tablename, fileName) {
    "use strict";
    $('#' + tablename).tableExport({
        type: 'json',
        fileName: fileName
    });
}

function exportListJSONTest() {
    "use strict";
    $("#offResultTable").tableExport({
        type: 'json'
    });
}

function stopCodon(cds, bestResult) {
    "use strict";
    $("#stopCodonBody").html("");
    let codonArr = cds.match(/\w{3}/ig);
    let stopArr = new Array();

    for (let i = 0; i < codonArr.length; i++) {
        if (codonArr[i].match(/(taa|tga|tag)/i)) {
            stopArr.push([i, codonArr[i]]);
        }
    }
    resultDis(stopArr, bestResult, cds.length);
    addStopCodonList(stopArr, cds.length);
}

function addStopCodonList(stopArr, len) {
    "use strict";
    let stopBody = document.getElementById("stopCodonBody");
    for (let i = 0; i < stopArr.length; i++) {
        let tr = document.createElement("tr");
        let td = document.createElement("td");
        td.innerHTML = i + 1;
        tr.appendChild(td);
        let td1 = document.createElement("td");
        let td2 = document.createElement("td");
        let td3 = document.createElement("td");
        td1.innerHTML = stopArr[i][1];
        td2.innerHTML = (stopArr[i][0] * 3 + 1);
        td3.innerHTML = Math.round((stopArr[i][0] * 3 + 1) / len * 10000) / 100 + "%";
        tr.appendChild(td1);
        tr.appendChild(td2);
        tr.appendChild(td3);
        stopBody.appendChild(tr);
    }
}



function showOrHide2() {
    "use strict";
    if (document.all.childAd2.style.display === 'none') {
        document.all.childAd2.style.display = "";
        document.getElementById("detailSetting2").innerHTML = "[-]Customize base editor settings:";
    } else {
        document.all.childAd2.style.display = "none";
        document.getElementById("detailSetting2").innerHTML = "[+]Customize base editor settings:";

    }
}

function generateBEselect(selectid, newselectid, tid) {
    $("#" + selectid).html("");
    let tmpsel = $('<select class="selectpicker form-control" id="' + newselectid + '"></select>');
    tmpsel.append($('<option>BE:"NGG,Spacer:20nt,Edit:[4-8]"</option>'));
    tmpsel.append($('<option>YE1-BE3:"NGG,Spacer:20nt,Edit:[5-7]"</option>'));
    tmpsel.append($('<option>EE-BE3:"NGG,Spacer:20nt,Edit:[5-6]"</option>'));
    tmpsel.append($('<option>YEE-BE3:"NGG,Spacer:20nt,Edit:[6-6]"</option>'));
    tmpsel.append($('<option>VQR-BE3:"NGAN,Spacer:20nt,Edit:[4-11]"</option>'));
    tmpsel.append($('<option>VRER-BE3:"NGCG,Spacer:20nt,Edit:[3-10]"</option>'));
    tmpsel.append($('<option>SaBE:"NNGRRT,Spacer:21nt,Edit:[3-12]"</option>'));
    tmpsel.append($('<option>Sa(KKH)-BE3:"NNNRRT,Spacer:21nt,Edit:[3-12]"</option>'));
    tmpsel.append($('<option>Cas12a–BE:"TTTV,Spacer:20nt,Edit:[10-12]"</option>'));
    tmpsel.append($('<option>Target-AID:"NGG,Spacer:20nt,Edit:[2-4]"</option>'));
    tmpsel.append($('<option>Target-AID-NG:"NG,Spacer:20nt,Edit:[2-4]"</option>'));
    tmpsel.append($('<option>xBE3:"NG,Spacer:20nt,Edit:[4-8]"</option>'));
    tmpsel.append($('<option>BE-PLUS:"NGG,Spacer:20nt,Edit:[4-14]"</option>'));
    $("#" + selectid).append(tmpsel);
    $("#" + newselectid).selectpicker();
    $("#" + newselectid).on("change", function () {
        selectModel(tid);
    });
}

function generateGeneticCodes(selectdiv, newselectid) {
    $("#" + selectdiv).html("");
    let tmpsel = $('<select class="selectpicker form-control" style="z-index:9999" id="' + newselectid + '"></select>');
    tmpsel.append($('<option>1. The Standard Code</option>'));
    tmpsel.append($('<option>2. The Vertebrate Mitochondrial Code</option>'));
    tmpsel.append($('<option>3. The Yeast Mitochondrial Code</option>'));
    tmpsel.append($('<option>4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code</option>'));
    tmpsel.append($('<option>5. The Invertebrate Mitochondrial Code</option>'));
    tmpsel.append($('<option>6. The Ciliate, Dasycladacean and Hexamita Nuclear Code</option>'));
    tmpsel.append($('<option>9. The Echinoderm and Flatworm Mitochondrial Code</option>'));
    tmpsel.append($('<option>10. The Euplotid Nuclear Code</option>'));
    tmpsel.append($('<option>11. The Bacterial, Archaeal and Plant Plastid Code</option>'));
    tmpsel.append($('<option>12. The Alternative Yeast Nuclear Code</option>'));
    tmpsel.append($('<option>13. The Ascidian Mitochondrial Code</option>'));
    tmpsel.append($('<option>14. The Alternative Flatworm Mitochondrial Code</option>'));
    tmpsel.append($('<option>16. Chlorophycean Mitochondrial Code</option>'));
    tmpsel.append($('<option>21. Trematode Mitochondrial Code</option>'));
    tmpsel.append($('<option>22. Scenedesmus obliquus Mitochondrial Code</option>'));
    tmpsel.append($('<option>23. Thraustochytrium Mitochondrial Code</option>'));
    tmpsel.append($('<option>24. Pterobranchia Mitochondrial Code</option>'));
    tmpsel.append($('<option>25. Candidate Division SR1 and Gracilibacteria Code</option>'));
    tmpsel.append($('<option>26. Pachysolen tannophilus Nuclear Code</option>'));
    tmpsel.append($('<option>27. Karyorelict Nuclear Code</option>'));
    tmpsel.append($('<option>28. Condylostoma Nuclear Code</option>'));
    tmpsel.append($('<option>29. Mesodinium Nuclear Code</option>'));
    tmpsel.append($('<option>30. Peritrich Nuclear Code</option>'));
    tmpsel.append($('<option>31. Blastocrithidia Nuclear Code</option>'));
    tmpsel.append($('<option>33. Cephalodiscidae Mitochondrial UAA-Tyr Code</option>'));
    $("#" + selectdiv).append(tmpsel);
    $("#" + newselectid).selectpicker();
    $("#" + newselectid).change(function () {
        setORFAtgTgaCodons(newselectid);
    });
    setORFAtgTgaCodons(newselectid);
}


function selectModel(num) {
    let mod = $("#PAMmodels" + num).val();
    let infoDraw = selectModelAndDraw(mod);
    $("#PAMValue" + num).val(infoDraw[1]);
    $("#spacerLenValue" + num).val(infoDraw[0] - infoDraw[1].length);
    $("#spacerLenSlider" + num).bootstrapSlider("setValue", infoDraw[0] - infoDraw[1].length);
    $("#spacerValue" + num).val(infoDraw[4] + "-" + infoDraw[5]);
    $("#spacerSlider" + num).bootstrapSlider("setValue", [infoDraw[4], infoDraw[5]]);
    if (mod === "Cas12a–BE") {
        $("#spacerLocate" + num).val("3' of PAM");
    } else {
        $("#spacerLocate" + num).val("5' of PAM");
    }
    drawPAM('figPAM' + num, infoDraw[0], infoDraw[1], infoDraw[2], infoDraw[3], [infoDraw[4], infoDraw[5]], infoDraw[6]);
}


function selectModelAndDraw(mod) {
    "use strict";
    let infoDraw = new Array(6);
    if (mod === 'BE:"NGG,Spacer:20nt,Edit:[4-8]"') {
        infoDraw = [23, "NGG", "NCC", 21, 4, 8, "+"];
    } else if (mod === 'YE1-BE3:"NGG,Spacer:20nt,Edit:[5-7]"') {
        infoDraw = [23, "NGG", "NCC", 21, 5, 7, "+"];
    } else if (mod === 'EE-BE3:"NGG,Spacer:20nt,Edit:[5-6]"') {
        infoDraw = [23, "NGG", "NCC", 21, 5, 6, "+"];
    } else if (mod === 'YEE-BE3:"NGG,Spacer:20nt,Edit:[6-6]"') {
        infoDraw = [23, "NGG", "NCC", 21, 6, 6, "+"];
    } else if (mod === 'VQR-BE3:"NGAN,Spacer:20nt,Edit:[4-11]"') {
        infoDraw = [24, "NGAN", "NCTN", 21, 4, 11, "+"];
    } else if (mod === 'VRER-BE3:"NGCG,Spacer:20nt,Edit:[3-10]"') {
        infoDraw = [24, "NGCG", "NCGC", 21, 3, 10, "+"];
    } else if (mod === 'SaBE:"NNGRRT,Spacer:21nt,Edit:[3-12]"') {
        infoDraw = [27, "NNGRRT", "NNCYYA", 22, 3, 12, "+"];
    } else if (mod === 'Sa(KKH)-BE3:"NNNRRT,Spacer:21nt,Edit:[3-12]"') {
        infoDraw = [27, "NNNRRT", "NNNYYA", 22, 3, 12, "+"];
    } else if (mod === 'Cas12a–BE:"TTTV,Spacer:20nt,Edit:[10-12]"') {
        infoDraw = [27, "TTTV", "CCCB", 1, 10, 12, "-"];
    } else if (mod === 'Target-AID:"NGG,Spacer:20nt,Edit:[2-4]"') {
        infoDraw = [23, "NGG", "NCC", 21, 2, 4, "+"];
    } else if (mod === 'Target-AID-NG:"NG,Spacer:20nt,Edit:[2-4]"') {
        infoDraw = [22, "NG", "NC", 21, 2, 4, "+"];
    } else if (mod === 'xBE3:"NG,Spacer:20nt,Edit:[4-8]"') {
        infoDraw = [22, "NG", "NC", 21, 4, 8, "+"];
    } else if (mod === 'BE-PLUS:"NGG,Spacer:20nt,Edit:[4-14]"') {
        infoDraw = [23, "NGG", "NCC", 21, 4, 14, "+"];
    }
    return infoDraw;

}

function demoSpacerSeq() {
    $("#spacerSequence2").val("TTACGCCAGCTGGCGAAAGG\nTATTACGCCAGCTGGCGAAA\nCAACAGTTGCGCAGCCTGAA\nCAAAGCGCCATTCGCCATTC\nCAGACGCGAATTATTTTTGA\nTTTATGGCAGGGTGAAACGC\nGGATGAGCAGACGATGGTGC\nGTGTACCACAGCGGATGGTT\nGCGACCAGATGATCACACTC\nCCGGTGCAGTATGAAGGCGG");
}

function getValueForDraw() {
    "use strict";
    let infoDraw = new Array(6);
    if ($("input[name='besetting']:checked").val() == 'sel') {
        infoDraw = selectModelAndDraw($("#PAMmodels1").val());
    } else {
        if ($("input[name='cuspam']:checked").val() == "sel") {
            infoDraw[1] = $("#PAMValue1").val();
        } else {
            infoDraw[1] = $("#cunpamtext").val();
        }
        infoDraw[2] = transStrand(infoDraw[1]);
        infoDraw[2] = infoDraw[2].split("").reverse().join("");
        infoDraw[0] = Number($("#spacerLenValue1").val()) + infoDraw[1].length;

        let spRegion = $("#spacerValue1").val().split("-");
        infoDraw[4] = Number(spRegion[0]);
        infoDraw[5] = Number(spRegion[1]);
        if ($("#spacerLocate1").val() === "3' of PAM") {
            infoDraw[6] = "-";
            infoDraw[3] = 1;
        } else {
            infoDraw[6] = "+";
            infoDraw[3] = infoDraw[0] - infoDraw[1].length + 1;
        }
    }
    return infoDraw;
}

function getValueForDraw2() {
    "use strict";
    let infoDraw = new Array(6);
    if ($("input[name='besetting2']:checked").val() == 'sel') {
        infoDraw = selectModelAndDraw($("#PAMmodels2").val());
    } else {
        if ($("input[name='cuspam2']:checked").val() == "sel2") {
            infoDraw[1] = $("#PAMValue2").val();
        } else {
            infoDraw[1] = $("#cunpamtext2").val().toUpperCase();
        }
        infoDraw[2] = transStrand(infoDraw[1]);
        infoDraw[2] = infoDraw[2].split("").reverse().join("");
        infoDraw[0] = Number($("#spacerLenValue2").val()) + infoDraw[1].length;

        let spRegion = $("#spacerValue2").val().split("-");
        infoDraw[4] = Number(spRegion[0]);
        infoDraw[5] = Number(spRegion[1]);
        if ($("#spacerLocate2").val() === "3' of PAM") {
            infoDraw[6] = "-";
            infoDraw[3] = 1;
        } else {
            infoDraw[6] = "+";
            infoDraw[3] = infoDraw[0] - infoDraw[1].length + 1;
        }
    }
    return infoDraw;
}

function selectModelsChange() {
    "use strict";
    let infoDraw = getValueForDraw();
    if (infoDraw[5] > infoDraw[0] - infoDraw[1].length) {
        $("#spacerSlider1").bootstrapSlider("setValue", infoDraw[0] - infoDraw[1].length);
        $("#spacerValue1").val(infoDraw[4] + "-" + infoDraw[5]);
    }
    drawPAM('figPAM1', infoDraw[0], infoDraw[1], infoDraw[2], infoDraw[3], [infoDraw[4], infoDraw[5]], infoDraw[6]);
}

function selectModelsChange2() {
    "use strict";
    let infoDraw = getValueForDraw2();
    if (infoDraw[5] > infoDraw[0] - infoDraw[1].length) {
        $("#spacerSlider2").bootstrapSlider("setValue", infoDraw[0] - infoDraw[1].length);
        $("#spacerValue2").val(infoDraw[4] + "-" + infoDraw[5]);
    }
    drawPAM('figPAM2', infoDraw[0], infoDraw[1], infoDraw[2], infoDraw[3], [infoDraw[4], infoDraw[5]], infoDraw[6]);
}

function demoSeq() {
    "use strict";
    /*jshint multistr: true */
    $("#cdsSequence").val(">lacZ\nATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAACTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGACGGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAATGGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCATCATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATGAAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTACACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGCATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGCGTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGGAATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGATCCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATTTGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCCATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCCCACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTATCCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGATGAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGCCAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAAGCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACCAGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGATGGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTGATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGCGTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAGTGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCATCTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAACCGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCGCTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACCCGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAAGCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCTCACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGATGGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCGGCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGATTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGATCTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGCGGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATCAGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAAGAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGGAGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAAAAATAA");
}

function offPredictReset2() {
    $("#spacerSequence2").val("");
    $("#PAMmodels2").val("BE1");
    $("#fastaUpload2").val("");
    $("#childAd2").hide();
    $("input[value='sel']").prop('checked', true);
    $("#offListBody2").html("");
    let name = 2;
    hideSomething("warnWrongZip" + name);
    hideSomething("warnWrongFile" + name);
    hideSomething("warnText" + name);
    hideSomething("warnFasta" + name);
    hideSomething("warnNotSupport" + name);
}

function offPredictReset() {
    $("#fastaUpload").val("");
    $("#offListBody").html("");
}

function drawTwoPie(allp, cbeip, allm, cbeim) {
    let staticpie = echarts.init(document.getElementById("twoPieStatic"));
    let staticpieoption = {
        toolbox: {
            left: '95%',
            top: 'center',
            orient: 'vertical',
            itemSize: 20,
            itemGap: 20,
            iconStyle: {
                normal: {
                    textPosition: 'left'
                },
                emphasis: {
                    textPosition: 'top'
                }
            },
            feature: {
                dataView: {
                    title: "Data",
                    readOnly: true,
                    lang: ['View', 'Close', 'Refresh']
                },
                saveAsImage: {
                    pixelRatio: 2,
                    name: "PAMStatic",
                    title: "Save"
                }
            }
        },
        tooltip: {
            trigger: 'item',
            formatter: "{a} <br/>{b}: {c} ({d}%)"
        },
        legend: {
            data: ['Editable sites (Plus)', 'CBEI sites (Plus)', 'Editable sites (Minus)', 'CBEI sites (Minus)']
        },
        series: [{
            name: 'Plus strand',
            type: 'pie',
            selectedMode: 'single',
            center: ['30%', '60%'],
            data: [{
                    value: allp,
                    name: "Editable sites (Plus)"
                },
                {
                    value: cbeip,
                    name: "CBEI sites (Plus)",
                    selected: true
                }
            ]
        }, {
            type: 'pie',
            name: 'Minus strand',
            selectedMode: 'single',
            center: ['75%', '60%'],
            data: [{
                    value: allm,
                    name: "Editable sites (Minus)"
                },
                {
                    value: cbeim,
                    name: "CBEI sites (Minus)",
                    selected: true
                }
            ]
        }]
    };
    staticpie.setOption(staticpieoption, true);
};


function drawPAMStaBar(data) {
    let rpkmchart = echarts.init(document.getElementById('posBarStatic'));
    let rpkmoption = {
        legend: {
            data: data.name,
            align: 'left'
        },
        toolbox: {
            left: '95%',
            top: 'center',
            orient: 'vertical',
            itemSize: 20,
            itemGap: 20,
            iconStyle: {
                normal: {
                    textPosition: 'left'
                },
                emphasis: {
                    textPosition: 'top'
                }
            },
            feature: {
                magicType: {
                    type: ['stack', 'tiled'],
                    title: {
                        stack: 'Stack',
                        tiled: 'Tiled'
                    }
                },
                dataView: {
                    title: "Data",
                    readOnly: true,
                    lang: ['View', 'Close', 'Refresh']
                },
                saveAsImage: {
                    pixelRatio: 2,
                    name: "CBEIPositionBarplot",
                    title: "Save"
                }
            }

        },
        dataZoom: [{
            type: 'inside'
        }, {
            type: 'slider'
        }, {
            show: true,
            yAxisIndex: 0,
            filterMode: 'empty',
            width: 30,
            height: '64%',
            showDataShadow: false,
            left: '91%'
        }],
        tooltip: {},
        xAxis: {
            data: xlabpos,
            silent: false,
            name: "ORF relative position (%)",
            nameLocation: "center",
            nameGap: 30,
            splitLine: {
                show: false
            }
        },
        yAxis: [{
            name: "Counts",
            nameLocation: "center",
            nameGap: 30
        }],
        series: [],
        animationEasing: 'elasticOut',
        animationDelayUpdate: function (idx) {
            return idx * 5;
        }
    };
    let datasets = new Array();
    for (let i = 0; i < data.name.length; i++) {
        datasets.push({
            name: data.name[i],
            type: 'bar',
            data: data.count[i],
            stack: 'dms',
            animationDelay: function (idx) {
                return idx * 10 + i * 10;
            }
        });
    };
    rpkmoption.series = datasets;
    rpkmchart.setOption(rpkmoption, true);
};