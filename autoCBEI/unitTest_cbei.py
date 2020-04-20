import unittest, re
from cbei import cbei


class TestCBEI(unittest.TestCase):

    beinfos = {
        "BE": ["NGG", 20, 4, 8, 5],
        "YE1-BE3": ["NGG", 20, 5, 7, 5],
        "EE-BE3": ["NGG", 20, 5, 6, 5],
        "YEE-BE3": ["NGG", 20, 6, 6, 5],
        "VQR-BE3": ["NGAN", 20, 4, 11, 5],
        "VRER-BE3": ["NGCG", 20, 3, 10, 5],
        "SaBE": ["NNGRRT", 21, 3, 12, 5],
        "Sa(KKH)-BE3": ["NNNRRT", 21, 3, 12, 5],
        "Cas12a–BE": ["TTTV", 20, 10, 12, 3],
        "Target-AID": ["NGG", 20, 2, 4, 5],
        "Target-AID-NG": ["NG", 20, 2, 4, 5],
        "xBE3": ["NG", 20, 4, 8, 5],
        "BE-PLUS": ["NGG", 20, 4, 14, 5]
    }

    def setUp(self):
        self.maxDiff = "None"

    def test_regexPlus(self):
        # The program first converts the PAM sequence into a
        # regular expression and then matches it in the target sequence.
        # The PAM sequence often contains degenerate bases.
        pam = "AWTRCSNGHBVD"
        tseq = "cgctgatttcgAWTRCSNGHBVDacacgtacgtatAATGCCCGWGRTcgtacgatccgagcgcgcgcgcgATTGCGKGYTARcagtcgatcgatctacgtcaATTGCGNGMCCK"
        tseqr = cbei.tranStrand(tseq)
        pamReg = cbei.getPAMReg(pam)
        pos = pamReg.finditer(tseq)
        poslist = []
        for tpos in pos:
            poslist.append(tpos.start())
        self.assertEqual(poslist, [11, 35, 70, 102])
    
    def test_addLabel(self):
        # raw Seq: AAAATTTTCC(c->t)CGGGGCCC(c->t)TTTTAAAA
        # labeled Seq: AAAA|{TTTTCC(c->t)C[GGGG]CCC(c->t)TTTT}|AAAA
        skipReg = re.compile(r',|\{|\}|\[|\]|\)|\-|\>|[c,t]')
        tseq="AAAATTTTCC(c->t)CGGGGCCC(c->t)TTTTAAAA"
        ttarget="AAAA|{TTTTCC(c->t)C[GGGG]CCC(c->t)TTTT}|AAAA"
        tseq=cbei.addLabel(tseq,12,"[",skipReg)
        tseq=cbei.addLabel(tseq,16,"]",skipReg)
        tseq=cbei.addLabel(tseq,24,"}|",skipReg)
        tseq=cbei.addLabel(tseq,4,"|{",skipReg) # add at pos right
        self.assertEqual(tseq,ttarget)

    def test_regexMinus(self):
        pam = "AWTRCSNGHBVD"
        tseq = "cgctgatttcgAWTRCSNGHBVDacacgtacgtatAATGCCCGWGRTcgtacgatccgagcgcgcgcgcgATTGCGKGYTARcagtcgatcgatctacgtcaATTGCGNGMCCKacgtacgatcg"
        tseqr = cbei.tranStrand(tseq)
        pamRegr = cbei.getPAMReg(cbei.tranStrand(pam))
        posr = pamRegr.finditer(tseqr)
        posrlist = []
        for tpos in posr:
            posrlist.append(tpos.start())
        self.assertEqual(posrlist, [11, 43, 78, 102])

    def test_cbeiPredictPlus1(self):
        # BE: NGG
        # GT{C,GT[T,TTA,(c->t)]aa,CGT,CGT,GAC,T}|gg,g
        tseq = "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTAcaaCGTCGTGACTgggAAAACCCTGGCGTTACCTAA"
        resDict = {'mode': 'Plus', 'pamSeq': 'ggg', 'spacerRange': '30-49', 'editRange': '33-37', 'pamRange': '50-52', 'pureSpacer': 'CGTTTTAcaaCGTCGTGACT',
                   'detailSpacer': 'GT{C,GT[T,TTA,(C->T)]AA,CGT,CGT,GAC,T}|GG,G', 'editpos': '37', 'location': 0.5138888888888888, 'editpatten': 'Ac'}
        endReg = re.compile(r'cag|caa|cga', re.I)
        pamReg = cbei.getPAMReg("NGG")
        pos = pamReg.finditer(tseq)
        plusRes = cbei.calSapcer(tseq, pos, self.beinfos["BE"], endReg, 'Plus')
        self.assertEqual(len(plusRes), 1)
        self.assertDictEqual(plusRes[0], resDict)

    def test_cbeiPredictPlus2(self):
        # BE: NGG
        # {CAA,[(C->T)AG,TT]G,CGC,AGC,CTG,AA}|T,GG
        # {CAA,(C->T[)AG,TT]G,CGC,AGC,CTG,AAT}|,GG
        import re
        tseq = "GATCGCCCTTCCCAACAGTTGCGcagCCTGAAtggCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAA"
        resDict = {'mode': 'Plus', 'pamSeq': 'tgg', 'spacerRange': '13-32', 'editRange': '16-20', 'pamRange': '33-35', 'pureSpacer': 'CAACAGTTGCGcagCCTGAA',
                   'detailSpacer': '{CAA,[(C->T)AG,TT]G,CGC,AGC,CTG,AAT}|,GG', 'editpos': '16', 'location': 0.2318840579710145, 'editpatten': 'AC'}
        endReg = re.compile(r'cag|caa|cga', re.I)
        pamReg = cbei.getPAMReg("NGG")
        pos = pamReg.finditer(tseq)
        plusRes = cbei.calSapcer(tseq, pos, self.beinfos["BE"], endReg, 'Plus')
        self.assertEqual(len(plusRes), 1)
        self.assertDictEqual(plusRes[0], resDict)

    def test_cbeiPredictMinus1(self):
        # BE: NGG
        # Single cbei site
        # C,CC|{T,TTC,GCC,AGC,TG[G,CGT],AAT,A}GC
        # GC{T,AT[T,ACG,(C->T)]CA,GCT,GGC,GAA,A}|GG,G
        import re
        tseq = "GATCGCCCTTCCCAACAGTTGCGcagCCTGAAtggCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAA"
        resDict = {'mode': 'Minus', 'pamSeq': 'AGG', 'spacerRange': '49-30', 'editRange': '46-42', 'pamRange': '27-29', 'pureSpacer': 'CAAAGCGCCATTCGccaTTC',
                   'detailSpacer': 'GG{C,AA[A,GCG,(C->T)]CA,TTC,GCC,ATT,CAG}|,G', 'editpos': '42', 'location': 0.6086956521739131, 'editpatten': 'GC'}
        endRegM = re.compile(r'cca', re.I)
        pamReg = cbei.getPAMReg("NGG")
        posR = pamReg.finditer(cbei.tranStrand(tseq))
        minusRes = cbei.calSapcer(cbei.tranStrand(
            tseq), posR, self.beinfos["BE"], endRegM, 'Minus')
        self.assertEqual(len(minusRes), 1)
        self.assertDictEqual(minusRes[0], resDict)

    def test_cbeiPredictMinus2(self):
        # BE: NGG
        # Double cbei sites
        # cc,c|{CT,TTC,GCC,AGC,t[gg,CGT],AAT}
        # {ATT,[ACG,(c->t)(c->t)]A,GCT,GGC,GAA,AG}|g,gg
        tseq = "ATGAGCGTTACCAAACTTAATCGCCCTTGGCACATcccTTTCGCCAGCTCCCGTAATAGCGAAGAGGCCTGA"
        resDict = {'mode': 'Minus', 'pamSeq': 'TGG', 'spacerRange': '33-14', 'editRange': '30-26',
                   'pamRange': '11-13', 'pureSpacer': 'GTGCCAAGGGCGATTAAGTT',
                   'detailSpacer': '{GTG,[(C->T)(C->T)A,AG]G,GCG,ATT,AAG,TT}|T,GG',
                   'editpos': '30,29', 'location': 0.4166666666666667, 'editpatten': 'GC,CC'}
        endRegM = re.compile(r'cca', re.I)
        pamReg = cbei.getPAMReg("NGG")
        posR = pamReg.finditer(cbei.tranStrand(tseq))
        minusRes = cbei.calSapcer(cbei.tranStrand(
            tseq), posR, self.beinfos["BE"], endRegM, 'Minus')
        self.assertEqual(len(minusRes), 1)
        self.assertDictEqual(minusRes[0], resDict)

    def test_cbeiPredictMinusPolyC(self):
        # BE: NGG
        # Muti cbei sites
        # PolyC
        tseq = "ATGAGCGTTACCCAACTTAATCGCCTTGCAGCACATccccccccTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCTGA"
        tdetails=['GC{T,AT[T,ACG,(C->T)]CA,GCT,GGC,GAA,A}|GG,G', 
                    '{ATT,[ACG,(C->T)(C->T)]A,GCT,GGC,GAA,AG}|G,GG', 
                    'A{TT,A[CG,(C->T)(C->T)A,]GCT,GGC,GAA,AGG}|,GGG', 
                    'AT{T,AC[G,(C->T)(C->T)A,G]CT,GGC,GAA,AGG,G}|GG,G', 
                    '{ACG,[(C->T)(C->T)A,GC]T,GGC,GAA,AGG,GG}|G,GG', 
                    'A{CG,C[(C->T)A,GCT,]GGC,GAA,AGG,GGG}|,GGG']
        resDict = {'mode': 'Minus', 'pamSeq': 'AGG', 'spacerRange': '49-30', 'editRange': '46-42',
                   'pamRange': '27-29', 'pureSpacer': 'CAAAGCGCCATTCGccaTTC',
                   'detailSpacer': 'GG{C,AA[A,GCG,(C->T)]CA,TTC,GCC,ATT,CAG}|,G', 'editpos': '42', 'location': 0.6086956521739131, 'editpatten': 'GC'}
        endRegM = re.compile(r'cca', re.I)
        pamReg = cbei.getPAMReg("NGG")
        posR = pamReg.finditer(cbei.tranStrand(tseq))
        minusRes = cbei.calSapcer(cbei.tranStrand(
            tseq), posR, self.beinfos["BE"], endRegM, 'Minus')
        
        details=[]
        for tm in minusRes:
            details.append(tm['detailSpacer'])
        self.assertEqual(len(minusRes), 6)
        self.assertEqual(details,tdetails)

    def test_transStrand(self):
        trans = cbei.tranStrand("AWTRCSNG")
        self.assertEqual(trans, "CNSGYAWT")


if __name__ == '__main__':
    unittest.main()
