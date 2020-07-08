# CRISPR-CBEI

CRISPR-CBEI is dedicated to designing potential spacers for cytosine base editor mediated gene inactivation. The editing parameters of cytosine base editors are entirely customizable, including PAM, editing window, and search region. CRISPR-CBEI also integrated with visual ORF identification and interactive statistical charts. Furthermore, this designing tool also provides a pure front-end off-target prediction searching within size-unlimited user-defined genome file.

CrisprCBEI is available in three versions to accommodate three usage scenarios: online and command line. Both versions have the function to design cytosine base editor mediated gene inactivation (CBEI). In addition, the web-based versions (both online and offline) also support ORF detection and off-target prediction on the web page. For the ‘online version,’ users can easily access through the website (["https://taolab.nwsuaf.edu.cn/crisprcbei/"](https://taolab.nwsuaf.edu.cn/crisprcbei/) or ["https://atlasbioinfo.github.io/CRISPR-CBEI/"](https://atlasbioinfo.github.io/CRISPR-CBEI/) ). 

![Figure1](./img/Figure1.png)
<center> Figure 1. Schematic overview of CRISPR-CBEI. 

(A)The calculation process from user input to designed available spacers. (B) Schematic diagram of the Base editing. (C)The workflow of CRISPR-CBEI. The left part is the design process of cytosine base editor mediated gene inactivation, and the right part is front-end off-target prediction. (D) The composition diagram of base editors.</center>

# Usage

## HTML version

You can download the files directly from the directory and open index.html.

```
├── css
│   ├── bootstrap.min.css
│   ├── style.css
│   └── ...
├── js
│   ├── align.js
│   ├── ORFfind.min.js
│   └── ...
├── data
│   ├── Test.fa
│   └── TestCBEI.fa
├── fonts
│   ├── ...
├── img
│   ├── BE.png
│   ├── ...
├── index.html
```

## autocbei

In addition to calculating a large number of CDSs, we have introduced a command-line version of the tool: autocbei ["https://github.com/atlasbioinfo/CRISPR-CBEI/tree/master/autocbei"](https://github.com/atlasbioinfo/CRISPR-CBEI/tree/master/autocbei). 

The detailed introduction of autocbei please see:
https://github.com/atlasbioinfo/CRISPR-CBEI/tree/master/autocbei

Bug report please submit issues or contat atlasbioin4@gmail.com.