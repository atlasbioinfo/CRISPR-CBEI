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


function hideORF() {
    $("#childAdORF").hide();
};

function showORF() {
    $("#childAdORF").show();
};

function showBEset() {
    $("#childAd").show();
    $("#PAMmodels1").prop("disabled", true);
};

function hideBEset() {
    $("#childAd").hide();
    $("#PAMmodels1").prop("disabled", false);
};

function showBEset2() {
    $("#childAd2").show();
    $("#PAMmodels2").prop("disabled", true);
};

function hideBEset2() {
    $("#childAd2").hide();
    $("#PAMmodels2").prop("disabled", false);
};

function selectRatio(idname) {
    $("input[value='" + idname + "']").prop('checked', true);
}

$("#cbeiTab a").click(function (e) {
    e.preventDefault();
    $(this).tab('show');
});
var regionSlid = $("#searchRegionSlider").bootstrapSlider().on("change", regionChange);
var spacerSlid = $("#spacerSlider1").bootstrapSlider().change(spacerChange);
var spacerLenSlid = $("#spacerLenSlider1").bootstrapSlider().on("change", spacerLenChange);
var ORFLenSlid = $("#ORFminlength").bootstrapSlider().on("change", ORFLenChange);
var ORFviewerChart = echarts.init(document.getElementById('orfViewerDiv'));
var spacerSlid2 = $("#spacerSlider2").bootstrapSlider().change(spacerChange2);
var spacerLenSlid2 = $("#spacerLenSlider2").bootstrapSlider().on("change", spacerLenChange2);

generateBEselect("addBEselect_Sub", "PAMmodels1", 1);
selectModel(1);

$('#cbeiTab a[href="#submitCDS"]').on('shown.bs.tab', function (e) {
    if ($("#PAMmodels1").length == 0) {
        generateBEselect("addBEselect_Sub", "PAMmodels1", 1);
        selectModel(1);
    }
});

$('#cbeiTab a[href="#Off-Target"]').on('shown.bs.tab', function (e) {
    if ($("#PAMmodels2").length == 0) {
        generateBEselect("addBEselect_Off", "PAMmodels2", 2);
        selectModel(2);
    }
});

generateGeneticCodes('addGeneticCodes', 'selectGenetic');

$("#PAMValue1").change(selectModelsChange);
$("#cunpamtext").change(selectModelsChange);
$("#cunpamtext").bind('input propertychange', function () {
    selectModelsChange();
});
$("#spacerLenValue1").change(selectModelsChange);
$("#spacerLocate1").change(selectModelsChange);
$("#spacerValue1").change(selectModelsChange);

$("#PAMValue2").change(selectModelsChange2);
$("#cunpamtext2").change(selectModelsChange2);
$("#cunpamtext2").bind('input propertychange', function () {
    selectModelsChange2();
});
$("#spacerLenValue2").change(selectModelsChange2);
$("#spacerLocate2").change(selectModelsChange2);
$("#spacerValue2").change(selectModelsChange2);

$('[data-toggle="tooltip"]').tooltip();
$('.dropdown-toggle').dropdown();
$("#childAdORF").hide();
$("#childAd").hide();
$("#childAd2").hide();
$("#offResultDir").hide();


genORFtableHeader();
genResultTableHeader();
$("#misNum2").selectpicker();

$("#ORFidentificationDiv").hide();
$("#CBEIpredictionDIV").hide();
$("#preResultDiv").hide();
$("#calResultDev").hide();

$("#helpnavdiv").affix();
$("#helpnavdiv").on('affixed-top.bs.affix', function () {
    $(this).children().children().first().addClass('active');
});
$("#resOffPre").click(function () {
    let sres = $("#resultTable").bootstrapTable('getSelections');
    let tspacer = new Array();
    for (s in sres) {
        tspacer.push(sres[s].purespacer);
    }
    let ttspacer = tspacer.join("\n");
    $("#spacerSequence2").val(ttspacer);
    $('#cbeiTab a[href="#Off-Target"]').tab('show');
    let infoDraw = getValueForDraw();
    if (infoDraw[5] > infoDraw[0] - infoDraw[1].length) {
        $("#spacerSlider2").bootstrapSlider("setValue", infoDraw[0] - infoDraw[1].length);
        $("#spacerValue2").val(infoDraw[4] + "-" + infoDraw[5]);
    }
    drawPAM('figPAM2', infoDraw[0], infoDraw[1], infoDraw[2], infoDraw[3], [infoDraw[4], infoDraw[5]],
        infoDraw[6]);
});
Back2top.init();
console.log("CRISPR-CBEI is " + document.readyState);