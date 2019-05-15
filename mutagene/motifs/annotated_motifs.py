
annotated_motifs = [
    {
        'name': 'APOBEC',
        'logo': 'T[C>K]W or W[G>M]A',
        'motif': 'TCW',
        'position': 1,
        'ref': 'C',
        'alt': 'K',
        'references': 'Starrett G. J. et al. Nature Communications (2016), Roberts S.A. et al. Nature Genetics (2013), Glaser A.P. et al. Oncotarget (2018)'
    },
    {
        'name': 'AGEING',
        'logo': '[C>T]G or C[G>A]',
        'motif': 'CG',
        'position': 0,
        'ref': 'C',
        'alt': 'T',
        'references': 'Jarvis M.C. et al. JNCL Natl Cancer Inst (2018), Ruzicka M. et al. Plos One (2017)'
    },
    {
        'name': 'UV Light',
        'logo': 'Y[C>T] or [G>A]R',
        'motif': 'YC',
        'position': 1,
        'ref': 'C',
        'alt': 'T',
        'references': 'Jarvis M.C. et al. JNCL Natl Cancer Inst (2018)'
    },
    {
        'name': 'AID',
        'logo': 'N[T>S]W or W[A>W]N',
        'motif': 'NTW',
        'position': 1,
        'ref': 'T',
        'alt': 'S',
        'references': 'Tsukamoto et al. Scientific Reports (2017)'
    },
    {
        'name': 'Tobacco',
        'logo': '[C>R]G',
        'motif': 'CG',
        'position': 0,
        'ref': 'C',
        'alt': 'R',
        'references': 'Alexandrov L.B. et al. Nature (2013)'
    },
    {
        'name': 'Aristolochic Acid (AA)',
        'logo': 'C[T>A]G',
        'motif': 'CTG',
        'position': 1,
        'ref': 'T',
        'alt': 'A',
        'references': 'Hoang M. L. et al. Science Translational Medicine (2013)'
    }
]

"""
Motif & References
APOBEC
tCw --> tTw or wGa --> wAa
* Starrett G. J. et al. Nature Communications (2016)
* Roberts S.A. et al. Nature Genetics (2013)
* Glaser A.P. et al. Oncotarget (2018)
tCw --> tGw or wGa --> wCa

AGEING
Cg --> Tg  or cG --> cA
tCg --> tTg
tCt --> tAt
* Jarvis M.C. et al. JNCL Natl Cancer Inst (2018)
* Ruzicka M. et al. Plos One (2017)

UV Light
(c/t)C --> (c/t)T or G(g/a)  --> A(g/a)
* Jarvis M.C. et al. JNCL Natl Cancer Inst (2018)

AID (activation-induced cytidine deaminase)
wrC --> wrG or Gyn --> Cyn
C-AID
wrCy --> wrTy/wrGy
* Rogozin et al. Scientific Reports (2016)
* Tsukamoto et al. Scientific Reports (2017)

Polymerase-N
nTw --> nGw or wAn --> wCn
nTw --> nCw or wAn --> wGn
* Tsukamoto et al. Scientific Reports (2017)

Chemotherapy-induced (e.g. cisplatin)
cT --> gA or Ag --> Tc
Tt --> Ct or aA --> aG
tGg --> tTg  or cCa --> cAa
cGg --> cTg or cCg --> cAg
* Liu D. et al. Nature Communications (2017)
* Findlay J.M. Nature Communications (2016)
* Huang K.K. et al. Scientific Reports (2016)

Aristolochic Acid (AA) Ð causes urothelial carcinoma
cTg --> cAg or cAg --> cTg
* Hoang M. L. et al. Science Translational Medicine (2013)

Tobacco Smoke
aGa --> aAa
tGc --> tTc
* Alexandrov L.B. et al. Nature (2013)
* Phillips D.H. Carcinogenesis (2002)
* Pfeifer G.P. et al Oncogene (2002)
* Pleasance E.D. et al. Nature (2010)
* Imielinski M. et al. Cell (2010)

Temozolomide
C(c/t) --> T(c/t)
* Alexandrov L.B. et al. Nature (2013)
* Johnson B.E. et al. Science (2014)
* Yip S. et al Clin Cancer res (2009)

Mismatch repair deficiency (MSI)
Cg --> Tg
Ct --> At
* Alexandrov L.B. et al. Nature (2013)

Overlap of AID and CpG
WRCG/CGYW with C --> G
* Rogozin et al. Scientific Reports (2016)
"""
