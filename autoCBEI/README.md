# autoCBEI.py

We provided command-line versions that can be used to design the cytosine base editor mediated gene inactivation for large amounts of CDS. 

We developed the 'autoCBEI' to automate the calculation of potential CBEI loci of the target CDSs and to perform statistics and calculations.

## 1. Install

AutoCBEI.py is based on Python scripts, and the following are the required dependencies.

```
	Python 3.X
	biopython
	matlibplot
	numpy
```
Once python3 is installed, you can install the required dependencies through PIP. In general, ‘matlibplot’ and ‘numpy’ are preloaded and do not need to be installed.

```
	pip install biopython matlibplot numpy
```

## 2 Usage

First copy the following files locally:
```bash
	cbei/
		__init__.py
		cbei.py
		stat.py
	autoCBEI.py
	#And CDSs file in fasta format for the required computation
	#Example (Download from Ensembl): 
	Bacillus_subtilis.ASM69118v1.cds.all.fa
```

Simply run the following command
```Python
python autoCBEI.py Bacillus_subtilis.ASM69118v1.cds.all.fa
```


## 3. Settings

You can customize the Base editor in the autoCBEI.py file:
```python
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
	# Or, add your own Base editors
    # "BE name" :[PAM, SpLength, EditBeg, EditEnd, Direction]
}
```

![BE](./BE.png)


> **PAM**: PAM sequence.

> **SpLength**: Spacer length.

> **EditBeg**: Edit windows begin.

> **EditEnd**: Edit windows end.

> **Direction**: 5 or 3. Spacer is at the 5 'end or 3' end of the PAM sequence. The example in the figure is the 5.




## 4 Output

### 4.1 Output message

The output message should be:
```
Input file: .\Bacillus_subtilis.ASM69118v1.cds.all.fa
Base editors:
        BE      PAM     Spacer  EditBegin       EditEnd Direction
        BE      NGG     20      4       8       5
        YE1-BE3 NGG     20      5       7       5
        EE-BE3  NGG     20      5       6       5
        YEE-BE3 NGG     20      6       6       5
        VQR-BE3 NGAN    20      4       11      5
        VRER-BE3        NGCG    20      3       10      5
        SaBE    NNGRRT  21      3       12      5
        Sa(KKH)-BE3     NNNRRT  21      3       12      5
        Cas12a–BE      TTTV    20      10      12      3
        Target-AID      NGG     20      2       4       5
        Target-AID-NG   NG      20      2       4       5
        xBE3    NG      20      4       8       5
        BE-PLUS NGG     20      4       14      5
Start calculating: BE
CBEI calculation of base editor "BE" , done!
Start calculating: YE1-BE3
CBEI calculation of base editor "YE1-BE3" , done!
Start calculating: EE-BE3
CBEI calculation of base editor "EE-BE3" , done!
Start calculating: YEE-BE3
CBEI calculation of base editor "YEE-BE3" , done!
Start calculating: VQR-BE3
CBEI calculation of base editor "VQR-BE3" , done!
Start calculating: VRER-BE3
CBEI calculation of base editor "VRER-BE3" , done!
Start calculating: SaBE
CBEI calculation of base editor "SaBE" , done!
Start calculating: Sa(KKH)-BE3
CBEI calculation of base editor "Sa(KKH)-BE3" , done!
Start calculating: Cas12a–BE
CBEI calculation of base editor "Cas12a–BE" , done!
Start calculating: Target-AID
CBEI calculation of base editor "Target-AID" , done!
Start calculating: Target-AID-NG
CBEI calculation of base editor "Target-AID-NG" , done!
Start calculating: xBE3
CBEI calculation of base editor "xBE3" , done!
Start calculating: BE-PLUS
CBEI calculation of base editor "BE-PLUS" , done!
Calculate complete!
Begin statistics...
####################
Pie charts for different BEs have been generated.
Path:
        ./CBEIPlot\[BE names]_Bacillus_subtilis.statPie.tiff
Transcript statistics and mapping completed.
Path:
        ./CBEIPlot/Bacillus_subtilis_CDSlength.tiff
        ./CBEIPlot/Bacillus_subtilis_GCstats.tiff
        ./CBEIPlot/Bacillus_subtilis_CodonUsage.tiff
        ./CBEIPlot/Bacillus_subtilis_CDSlength.tiff
The comparison of CBEI ratio of different BE has been completed.
Path:
        ./CBEIPlot/Bacillus_subtilis.statBar.tiff
        ./CBEIPlot/Bacillus_subtilis.statROC.tiff
CBEI statistics complete
```
### 4.2 Output files

The output directory or files shold be:
```
	CBEIRaw/
		BE_Bacillus_subtilis.ASM69118v1.cds.all.fa.cbei
        xxxx.cbei
    CBEIRes/
        BE_Bacillus_subtilis.Thre05.tsv
        xxxx.Threxx.tsv
    CBEIPlot/
        xxxx.tiff
```

#### 4.2.1 .cbei files format

```
KDE22635	Minus	0.8674698795180723	{TGC,[TTT,C(C->T)]A,TAA,GAT,TAA,AA}|T,G	TGCTTTCCATAAGATTAAAA	222-203	216,215	219-215	TG	201-202	TC,CC
```

>1. Fasta title: E.g., lacZ
>2. Strand: E.g., Plus
>3. Relative position: Edit position/Gene length. E.g., 0.012032520325203253
>4. Rich info: E.g., GT{C,GT[T,TTA,(C->T)]AA,CGT,CGT,GAC,T}|GG,G
>5. Spacer: E.g., CGTTTTACAACGTCGTGACT
>6. Spacer position: E.g., 30-49
>7. Edit position: E.g., 37
>8. Edit windows: E.g., 33-37
>9. PAM:E.g., GGG
>10. PAM position: E.g., 50-52
>11. Edit pattern: The edit pattern indicated the adjacent nucleotide at 5’ of the editable cytosine. Typically, the in vitro activity of the base editors follows TC ≥ CC ≥ AC > GC. E.g., AC

#### 4.2.2 .tsv files with threshold

Filter according to different thresholds. 
>* Thre25: CBEI sites at the upsteam 25% of the CDS.
>* Thre50: CBEI sites at the upsteam 50% of the CDS.
>* Thre75: CBEI sites at the upsteam 75% of the CDS.

```
#Editable sites in the first 50% of the transcript
#Transcript	Strand	Position	CBEIdetail	Spacer	SpacerRegion	EditPosition	EditWindowsRegion	PAM	PAMregion	Pattern

KDE22636	Plus	0.3686274509803922	GG{A,TT[T,CTT,(C->T)]AA,CAA,ACA,GTA,C}|TG,G	ATTTCTTCAACAAACAGTAC	87-106	94	90-94	TGG	107-109	TC
```

### 4.3 Figures

```
Pie charts for different BEs have been generated.
Path:
        ./CBEIPlot/[BE names]_Bacillus_subtilis.statPie.tiff
Transcript statistics and mapping completed.
Path:
        ./CBEIPlot/Bacillus_subtilis_CDSlength.tiff
        ./CBEIPlot/Bacillus_subtilis_GCstats.tiff
        ./CBEIPlot/Bacillus_subtilis_CodonUsage.tiff
        ./CBEIPlot/Bacillus_subtilis_CDSlength.tiff
The comparison of CBEI ratio of different BE has been completed.
Path:
        ./CBEIPlot/Bacillus_subtilis.statBar.tiff
        ./CBEIPlot/Bacillus_subtilis.statROC.tiff
CBEI statistics complete
```
