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

var labelStr = new Array('+1', "+2", "+3", "-1", "-2", "-3");
var labeldict = {'+1':0, "+2":1, "+3":2, "-1":3, "-2":4, "-3":5};

function ORFtest() {
    
    $("#resultli").removeClass("disabled");
    $("#resultli a").attr("href", "#Result");
    $('#cbeiTab a[href="#Result"]').tab('show');
    
    let seqData = getCDS();
    for (tnameindex in seqData.seqName){
        if (seqData.seqName[tnameindex].substr(0,1) == ">"){
            seqData.seqName[tnameindex]=seqData.seqName[tnameindex].substr(1);
        }
        if (seqData.seqName[tnameindex].length>40){
            seqData.seqName[tnameindex]=seqData.seqName[tnameindex].substr(0,40);
        }
    }
    // seq,seqName，后面跟数组
    let augReg;
    let tgaReg;
    let mLength = 75;
    let intronArr=new Array();
    if ($("input[name='intron']:checked").val() == 'cusintron') {
        intronArr = getInputCDSRegion();
    }
    if ($("input[name='ORFfinder']:checked").val() == "ORFoption1") {
        augReg = new RegExp("A[T|U]G", "i");
        tgaReg = new RegExp("[T|U]AA|[T|U]AG|[T|U]GA", "i");
    } else {
        let tLength = Number($("#ORFLenValue").val());
        if (tLength === 0) {
            mLength = 75;
        } else {
            mLength = tLength;
        }
        let allReg=getAUGreg();
        augReg=allReg.augReg;
        tgaReg=allReg.tgaReg;
    }
    generatSelectList("orfselectdiv","orfselectlist",seqData,mLength, augReg, tgaReg, intronArr);
}

function generatSelectList(divname,newdivname,seqData,mLength, augReg, tgaReg, intronArr){
    $("#" + divname).html("");
    let tmpsel = $('<select class="selectpicker form-control"  data-live-search="true" id="' + newdivname + '"></select>');
    for (i in seqData.seqName){
        tmpsel.append($('<option>'+seqData.seqName[i]+'</option>')); 
    }
    $("#" + divname).append(tmpsel);
    $("#" + newdivname).selectpicker();
    $("#" + newdivname).on("change", function () {
        geneOneORF(newdivname,seqData,mLength, augReg, tgaReg, intronArr);
    });
    geneOneORF(newdivname,seqData,mLength, augReg, tgaReg, intronArr);
}

function geneOneORF(selectname,seqData,mLength, augReg, tgaReg, intronArr){
    let selectedFasta=$('#'+selectname).val();
    ORFfinderATG(seqData, mLength, augReg, tgaReg, intronArr,selectedFasta);
    $("#ORFidentificationDiv").show();
    $("#CBEIpredictionDIV").hide();
    $("#preResultDiv").hide();
    $("#calResultDev").hide();
    $("#helpnavdiv").affix();

}

function getAUGreg(){
    let augregarr=new Array();
    $("input[name='augcodons']:checked").each(function (){
       augregarr.push($(this).val());
    });
    let tgaregarr=new Array();
    $("input[name='tgacodons']:checked").each(function (){
        tgaregarr.push($(this).val());
     });
     return {'augReg': new RegExp(augregarr.join("|"),"i"),'tgaReg':new RegExp(tgaregarr.join("|"),"i")}
}

function genORFtableHeader() {
    $("#orfTable").bootstrapTable({
        striped: true,
        pagination: true,
        sidePagination: "client",
        pageNumber: 1,
        pageSize: 10,
        strictSearch: true,
        showColumns: true,
        clickToSelect: true,
        uniqueId: "id",
        showExport: true,
        exportDataType: "all",
        columns: [{
            field: 'id',
            title: 'ID',
            sortable: true,
            width: "80"
        }, {
            field: 'fname',
            title: 'Title',
            sortable: true
        }, {
            field: 'frame',
            title: 'Frame',
            sortable: true
        }, {
            field: 'start',
            title: 'Start',
            sortable: true
        }, {
            field: 'end',
            title: 'End',
            sortable: true
        }, {
            field: 'length',
            title: 'Length',
            sortable: true
        }, {
            field: 'orfSeq',
            title: 'Sequence',
            width: "250"
        }, {
            field: 'predict',
            title: 'Predict',
            align: 'center',
            // events: window.predictEvents,
            formatter: predictFormatter
        }],
        formatSearch: function () {
            return 'Search Fasta title';
        }
    });
}

function predictFormatter(value, row, index) {
    return '<a type="button" href="#CBEIpredictionDIV" class="btn btn-primary" onclick=\'showORFinfo("' + row.fname + '","' + row.frame + '",' + row.start + ',' + row.end + ',"' + row.orfSeq + '")\'>Predict</a>';
}

function showORFinfo(title, frame, start, end, orfSeq) {

    $("#CBEIpredictionDIV").show();
    location.href="#CBEIpredictionDIV";
    $("#preResultDiv").hide();
    $("#calResultDev").hide();

    let $table = $("#orfInfoTable");
    $table.bootstrapTable({
        striped: true,
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
        'item': 'Title',
        'value': title
    });
    subcontain.push({
        'item': 'Frame',
        'value': frame
    });
    subcontain.push({
        'item': 'Start',
        'value': start
    });
    subcontain.push({
        'item': 'End',
        'value': end
    });
    subcontain.push({
        'item': 'Length',
        'value': (end - start + 1)
    });
    subcontain.push({
        'item': 'Codons',
        'value': Math.abs((end - start + 1) / 3)
    });
    subcontain.push({
        'item': 'Sequence',
        'value': orfSeq
    });
    $table.bootstrapTable('load', subcontain);
    // calSpacer(title,frame,start,end,orfSeq);
}


function ORFfinderATG(seqData, mLength, augReg,tgaReg,intronArr,selectedFasta) {
    let tORFres = new Array();
    let index = 1;
    let maxL=-1;
    let target=0;
    for (let i=0;i<seqData.seqName.length;i++){
        if (seqData.seqName[i] == selectedFasta){
            target=i;
            break;
        }
    }
    
        let fname = seqData.seqName[target];
        if (!fname) {
            fname = "Undefined";
        }
        let seq = seqData.seq[target];
        if (intronArr.length != 0){
            seq=getIntron(seq, intronArr);
        }
        if (seq.length>maxL){
            maxL=seq.length;
        }
        let readingFrames = [1, 2, 3, -1, -2, -3];
        let rseq = transStrand(seq);
        for (readingFrame in readingFrames) {
            if (readingFrames[readingFrame] < 0) {
                tseq = rseq;
            } else {
                tseq = seq;
            }
            for (let i = Math.abs(readingFrames[readingFrame]) - 1; i < tseq.length - 2; i += 3) {
                if (augReg.test(tseq.substring(i, i + 3))) {
                    let lastCodon = tseq.length - (tseq.length % 3);
                    for (let j = i + 3; j <= lastCodon; j += 3) {
                        let tcodon = tseq.substring(j, j + 3);
                        if (tgaReg.test(tcodon) || j === lastCodon) {
                            if (j + 3 - i >= mLength) {
                                tORFres.push(ORFout(index, i, j, j === lastCodon, fname, readingFrames[readingFrame], tseq));
                                index++;
                            }
                            i = j;
                            break;
                        }
                    }
                }
            }
        }
    
    index--;
    if (seqData.seqName.length == 0) {
        $("#ORFtoolbarText").html('<strong class="ORFtoolbarTextstrong">' + index + '</strong> ORFs.');
    } else {
        $("#ORFtoolbarText").html('<strong class="ORFtoolbarTextstrong">' + seqData.seqName.length + '</strong> Fasta sequence, <strong class="ORFtoolbarTextstrong">' + index + '</strong> ORFs of <strong>\"'+selectedFasta+ "\"</strong>.");
    }
    if (tORFres.length == 0) {
        tORFres.push(ORFout(1, 0, seqData.seq[0].length - 3, true, seqData.seqName, 1, seqData.seq[0]));
        showORFWarn(2);
    } else {
        showORFWarn(1);
    }
    $("#orfTable").bootstrapTable("load", tORFres);
    let orfInfo=new Array();
    for (let i in tORFres){
        orfInfo.push({"name":tORFres[i].id,"label":tORFres[i].frame,"beg":tORFres[i].start,"end":tORFres[i].end,'fasta':tORFres[i].fname,'seq':tORFres[i].orfSeq});
    }
    
    drawORFviewer(orfInfo,maxL,selectedFasta);
}

function getIntron(seq, intronArr){
    let newSeq="";
    for (let i=0;i<intronArr.length;i+=2){
        newSeq=newSeq+seq.substring(intronArr[i]-1,intronArr[i+1]);
    }
    return newSeq;
}

function showORFWarn(num) {
    // class="alert alert-danger"
    // style="visibility: hidden;display: none"
    if (num == 1) {
        $("#warnORFdiv").show();
        $("#warnORFdiv").html('<div class="alert alert-success fade in"><button class="close" type="button" data-dismiss="alert">&times;</button>ORF identification success !</div>')
    } else if (num == 2) {
        $("#warnORFdiv").show();
        $("#warnORFdiv").html('<div class="alert alert-danger fade in"><button class="close" type="button" data-dismiss="alert">&times;</button>Unrecognized ORF. If you are sure you need to predict non-ORF, please continue!</div>')
    }
}

function ORFout(index, i, j, lastCodon, fname, rframe, seq) {
    let start = i + 1;
    let stop = j + 3;
    let length = j + 3 - i;
    let sequence = seq.substring(i, j + 3);
    if (rframe < 0) {
        start = seq.length - start + 1;
        stop = seq.length - stop + 1;
    }
    // if (lastCodon === true) {
    //     stop = ">" + stop;
    //     length = length + "+"
    // }
    if (rframe > 0) {
        rframe = "+" + rframe;
    }
    let tbutt = '<button class="btn btn-primary">Predict</button>';
    return {
        "id": "ORF"+index,
        "fname": fname,
        "frame": rframe,
        "start": start,
        "end": stop,
        "length": length,
        "orfSeq": sequence,
        "Predict": tbutt
    };
}

function customSearch(data, text) {
    return data.filter(function (row) {
        return row.fname.indexOf(text) > -1
    })
}

function setORFAtgTgaCodons(newid) {

    $("#ORFidentificationDiv").hide();

    let selectVal = $("#"+newid).val();
    $("#orfStartCodons").html("");
    $("#orfStopCodons").html("");
    let augcodon = new Array();
    let tgacodon = new Array();

    switch (selectVal) {
        case '1. The Standard Code':
            augcodon = ['ATG', "GTG", "TTG", "CTG"];
            tgacodon = ['TAA', 'TAG', 'TGA'];
            break;
        case '12. The Alternative Yeast Nuclear Code':
        case '26. Pachysolen tannophilus Nuclear Code':
            augcodon = ['ATG', "CTG"];
            tgacodon = ['TAA', 'TAG', 'TGA'];
            break;
        case '22. Scenedesmus obliquus Mitochondrial Code':
        case '28. Condylostoma Nuclear Code':
            augcodon = ['ATG'];
            tgacodon = ['TAA', 'TAG', 'TGA'];
            break;
        case '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code':
        case '27. Karyorelict Nuclear Code':
        case '29. Mesodinium Nuclear Code':
        case '30. Peritrich Nuclear Code':
            augcodon = ['ATG'];
            tgacodon = ['TGA'];
            break;
        case '2. The Vertebrate Mitochondrial Code':
            augcodon = ['ATG', 'ATA', 'ATC', 'ATT', 'GTG'];
            tgacodon = ['TAA', 'TAG', 'AGA', 'AGG'];
            break;
        case '3. The Yeast Mitochondrial Code':
            augcodon = ['ATG', 'ATA', 'GTG'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code':
            augcodon = ['ATG', 'ATA', 'ATC', 'ATT', 'GTG', 'CTG', 'TTG', 'TTA'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '5. The Invertebrate Mitochondrial Code':
            augcodon = ['ATG', 'ATA', 'ATC', 'ATT', 'GTG', 'TTG'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '9. The Echinoderm and Flatworm Mitochondrial Code</option>':
            augcodon = ['ATG', 'GTG'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '10. The Euplotid Nuclear Code':
            augcodon = ['ATG'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '11. The Bacterial, Archaeal and Plant Plastid Code':
            augcodon = ['ATG', 'ATA', 'ATC', 'ATT', 'GTG', 'CTG', 'TTG'];
            tgacodon = ['TAA', 'TAG', 'TGA'];
            break;
        case '13. The Ascidian Mitochondrial Code':
            augcodon = ['ATG', "ATA", "TTG", "GTG"];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '14. The Alternative Flatworm Mitochondrial Code':
            augcodon = ['ATG'];
            tgacodon = ['TAG'];
            break;
        case '16. Chlorophycean Mitochondrial Code':
            augcodon = ['ATG'];
            tgacodon = ['TAA', 'TGA'];
            break;
        case '21. Trematode Mitochondrial Code':
            augcodon = ['ATG', "GTG"];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '23. Thraustochytrium Mitochondrial Code':
            augcodon = ['ATG', "ATT", "GTG"];
            tgacodon = ['TAA', 'TAG', 'TGA', 'TTA'];
            break;
        case '24. Pterobranchia Mitochondrial Code':
            augcodon = ['ATG', "GTG", "TTG", "CTG"];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '25. Candidate Division SR1 and Gracilibacteria Code':
            augcodon = ['ATG', "GTG", "TTG"];
            tgacodon = ['TAA', 'TAG'];
            break;

        case '31. Blastocrithidia Nuclear Code':
            augcodon = ['ATG'];
            tgacodon = ['TAA', 'TAG'];
            break;
        case '33. Cephalodiscidae Mitochondrial UAA-Tyr Code':
            augcodon = ['ATG', "GTG", "TTG", "CTG"];
            tgacodon = ['TAG'];
            break;
        default:
            augcodon = ['ATG', "GTG", "TTG", "CTG"];
            tgacodon = ['TAA', 'TAG', 'TGA'];
    }
    for (let startC in augcodon) {
        if (/A[T|U]G/.test(augcodon[startC])) {
            $("#orfStartCodons").append('<label class="checkbox-inline"><input type="checkbox" name="augcodons" value="' + augcodon[startC] + '" checked>' + augcodon[startC] + '</label>');
        } else {
            $("#orfStartCodons").append('<label class="checkbox-inline"><input type="checkbox" name="augcodons" value="' + augcodon[startC] + '">' + augcodon[startC] + '</label>');
        }

    }
    for (let stopC in tgacodon) {
        $("#orfStopCodons").append('<label class="checkbox-inline"><input type="checkbox" name="tgacodons" value="' +tgacodon[stopC] + '" checked>' + tgacodon[stopC] + '</label>');
    }
}

function drawORFviewer(orfInfo, geneLength, fastaName) {

    let data = [];
    let categories = ['+1', '+2', '+3', '-1', '-2', '-3'];
    let cColors = ['#d64f44', '#fcaf17', '#f173ac', '#007d65', '#009ad6', '#999d9c']

    for (let i in orfInfo) {
        let tlab=labeldict[orfInfo[i].label];
        if (orfInfo[i].beg>orfInfo[i].end){
            let tt=orfInfo[i].beg;
            orfInfo[i].beg=orfInfo[i].end;
            orfInfo[i].end=tt;
        }
        data.push({
            name: orfInfo[i].name,
            value: [
                tlab,
                orfInfo[i].beg,
                orfInfo[i].end,
                orfInfo[i].fasta,
                orfInfo[i].seq
            ],
            itemStyle: {
                normal: {
                    color: cColors[tlab]
                }
            }
        });
    }

    function renderItem(params, api) {

        let label = api.value(0);
        let beg = api.coord([api.value(1), label]);
        let end = api.coord([api.value(2), label]);
        let height = api.size([0, 1])[1] * 0.4;

        let rectShape = echarts.graphic.clipRectByRect({
            x: beg[0],
            y: beg[1] - height / 2,
            width: end[0] - beg[0],
            height: height
        }, {
            x: params.coordSys.x,
            y: params.coordSys.y,
            width: params.coordSys.width,
            height: params.coordSys.height
        });

        return rectShape && {
            type: 'rect',
            shape: rectShape,
            style: api.style()
        };
    }


    let option = {
        tooltip: {
            triggerOn: 'mousemove',
            formatter: function (params) {
                if (params.value[3].length>10){
                    params.value[3]=params.value[3].substr(0,10)+'...';
                }
                return [params.marker + params.name,
                    'Frame: ' + labelStr[params.value[0]],
                    'Fasta: ' + params.value[3],
                    'Range: ' + params.value[1] + 'nt - ' + params.value[2] + 'nt',
                    'Length: ' + (params.value[2] - params.value[1] + 1) + 'nt'
                ].join('</br>');
            }
        },
        title: {
            text: 'ORF viewer of ' + fastaName,
            left: 'center'
        },
        dataZoom: [{
            type: 'slider',
            filterMode: 'weakFilter',
            showDataShadow: false,
            top: 300,
            height: 10,
            borderColor: 'transparent',
            backgroundColor: '#e2e2e2',
            handleIcon: 'M10.7,11.9H9.3c-4.9,0.3-8.8,4.4-8.8,9.4c0,5,3.9,9.1,8.8,9.4h1.3c4.9-0.3,8.8-4.4,8.8-9.4C19.5,16.3,15.6,12.2,10.7,11.9z M13.3,24.4H6.7v-1.2h6.6z M13.3,22H6.7v-1.2h6.6z M13.3,19.6H6.7v-1.2h6.6z', // jshint ignore:line
            handleSize: 20,
            handleStyle: {
                shadowBlur: 6,
                shadowOffsetX: 1,
                shadowOffsetY: 2,
                shadowColor: '#aaa'
            },
            labelFormatter: ''
        }, {
            type: 'inside',
            filterMode: 'weakFilter'
        }],
        grid: {
            height: 200
        },
        xAxis: {
            min: 1,
            max: geneLength,
            scale: true,
            axisLabel: {
                formatter: function (val) {
                    return val + ' nt';
                }
            }
        },
        yAxis: {
            data: categories
        },
        series: [{
            type: 'custom',
            renderItem: renderItem,
            itemStyle: {
                normal: {
                    opacity: 0.8
                }
            },
            encode: {
                x: [1, 2],
                y: 0
            },
            data: data
        }]
    };

    ORFviewerChart.setOption(option, true);
    ORFviewerChart.on('click', function (params) {
        
        showORFinfo(params.data.value[3],labelStr[params.data.value[0]],params.data.value[1],params.data.value[2],params.data.value[4]);
    });
}