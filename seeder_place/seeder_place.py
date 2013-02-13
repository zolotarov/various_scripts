def seeder_result_analysis(seeder_result_file):
    """
    This script will search the file containing motifs found by Seeder and will create an
    HTML page for motifs that are found on PLACE
    
    example seeder_result_file content:
    10009 CGTGGC 3.62843e-03
    10015 CCGACC 7.29131e-03
    10033 AGGCGG 1.25627e-04
    10033 TATATA 2.51254e-04
    """
    import re
    from Bio.Seq import Seq

    # the dictionary contains the motif name as a key and the regular expression
    # for that motif and it's web address as values in a list
    PLACE_db = {"WINPSTPIIIK": ["AAGCGTAAGT", "S000495"],
                    "NODCON1GM": ["AAAGAT", "S000461"],
                    "AACACOREOSGLUB1": ["AACAAAC", "S000353"],
                    "BIHD1OS": ["TGTCA", "S000498"],
                    "SITEIIATCYTC": ["TGGGC[C|T]", "S000474"],
                    "GAREAT": ["TAACAA[A|G]", "S000439"],
                    "CAREOSREP1": ["CAACTC", "S000421"],
                    "14BPATERD1": ["CACTAAATTGTCAC", "S000412"],
                    "WBOXATNPR1": ["TTGAC", "S000390"],
                    "EREGCC": ["TAAGAGCCGCC", "S000036"],
                    "E2F1OSPCNA": ["GCGGGAAA", "S000396"],
                    "MYB1LEPR": ["GTTAGTT", "S000443"],
                    "PIATGAPB": ["GTGATCAC", "S000381"],
                    "ACIIPVPAL2": ["CCACCAACCCCC", "S000193"],
                    "SORLIP3AT": ["CTCAAGTGA", "S000484"],
                    "ABRELATERD1": ["ACGTG", "S000414"],
                    "CURECORECR": ["GTAC", "S000493"],
                    "CMSRE1IBSPOA": ["TGGACGG", "S000511"],
                    "MARCEN3": ["TGTTT[A|T]TG.TTTCCGAAA....[A|T][A|T][A|T]", "S000065"],
                    "CELLCYCLESC": ["CACGAAAA", "S000031"],
                    "MYCATRD22": ["CACATG", "S000174"],
                    "RSEPVGRP18": ["CATCCAACTTTCATATCCATGTGCTT", "S000289"],
                    "SORLIP5AT": ["GAGTGAG", "S000486"],
                    "HY5AT": ["TGACACGTGGCA", "S000345"],
                    "BOX2PSGS2": ["TCTAAGCAAAG", "S000204"],
                    "RGATAOS": ["CAGAAGATA", "S000191"],
                    "GCCCORE": ["GCCGCC", "S000430"],
                    "ATRICHPSPETE": ["AATATACTAGTATTATTTACTAAAAAAAATC", "S000248"],
                    "GCBP2ZMGAPC4": ["GTGGGCCCG", "S000351"],
                    "MSACRCYM": ["AGACCGTTG", "S000236"],
                    "MARARS": ["[A|T]TTTAT[A|G]TTT[A|T]", "S000064"],
                    "MREATCHS": ["TCTAACCTACCA", "S000356"],
                    "PE2FNTRNR1A": ["ATTCGCGC", "S000455"],
                    "BOX2PVCHS15": ["CTGTGATTAAATAT", "S000209"],
                    "CARG1ATAP3": ["GTTTACATAAATGGAAAA", "S000347"],
                    "HEXAT": ["TGACGTGG", "S000334"],
                    "MRNA3ENDTAH3": ["AATGGAAATG", "S000069"],
                    "WBOXNTCHN48": ["CTGAC[C|T]", "S000508"],
                    "CEREGLUBOX3PSLEGA": ["TGTAAAAGT", "S000034"],
                    "AT1BOX": ["AATATTTTTATT", "S000025"],
                    "L1DCPAL1": ["ATTCACCTACCC", "S000504"],
                    "SURECOREATSULTR11": ["GAGAC", "S000499"],
                    "POLLEN1LELAT52": ["AGAAA", "S000245"],
                    "RBENTGA3": ["TCCAACTTGGA", "S000358"],
                    "SORLREP4AT": ["CTCCTAATT", "S000489"],
                    "HSELIKENTACIDICPR1": ["C..GAA...TTC..G", "S000056"],
                    "C1MOTIFZMBZ2": ["TAACT[G|C]AGTTA", "S000237"],
                    "RYREPEATVFLEB4": ["CATGCATG", "S000102"],
                    "TAAAGSTKST1": ["TAAAG", "S000387"],
                    "SURE1STPAT21": ["AATAGAAAA", "S000186"],
                    "CCTCGTGTCTCGMGH3": ["CCTCGTGTCTC", "S000369"],
                    "ELEMENT1GMLBC3": ["GATATATTAATATTTTATTTTATA", "S000319"],
                    "AMMORESIIUDCRNIA1": ["GG[A|T]AGGGT", "S000374"],
                    "TATABOX1": ["CTATAAATAC", "S000108"],
                    "TACBBFNTEAS4": ["ACTCTACAGTACTC", "S000257"],
                    "AGMOTIFNTMYB2": ["AGATCCAA", "S000444"],
                    "S1FSORPL21": ["ATGGTATT", "S000215"],
                    "UPRMOTIFIIAT": ["CC............CCACG", "S000426"],
                    "HBOXPVCHS15": ["CCTACC.......CT....A", "S000361"],
                    "GRAZMRAB17": ["CACTGGCCGCCC", "S000150"],
                    "VOZATVPP": ["GCGT.......ACGC", "S000456"],
                    "DOFCOREZM": ["AAAG", "S000265"],
                    "SORLREP5AT": ["TTGCATGACT", "S000490"],
                    "ABREMOTIFIOSRAB16B": ["AGTACGTGGC", "S000290"],
                    "5256BOXLELAT5256": ["TGTGGTTATATA", "S000279"],
                    "SPHZMC1": ["CGTCCATGCAT", "S000293"],
                    "AGATCONSENSUS": ["TT[A|T]CC[A|T][A|T][A|T][A|T]..GG[A|T][A|T]", "S000316"],
                    "ACGTABREMOTIFA2OSEM": ["ACGTG[G|T]C", "S000394"],
                    "OCSGMHSP26A": ["TGATGTAAGAGATTACGTAA", "S000357"],
                    "TRANSTART": ["TAAACAATGGCT", "S000113"],
                    "GLUTEBOX1OSGT2": ["ATATCATGAGTCACTTCA", "S000046"],
                    "GLUTEBOX1OSGT3": ["TATCTAGTGAGTCACTTCA", "S000128"],
                    "E2FANTRNR": ["TTTCCCGC", "S000366"],
                    "PYRIMIDINEBOXHVEPB1": ["TTTTTTCC", "S000298"],
                    "CACGCAATGMGH3": ["CACGCAAT", "S000368"],
                    "GTGANTG10": ["GTGA", "S000378"],
                    "TATABOX2": ["TATAAAT", "S000109"],
                    "GAREHVAMY1": ["GGCCGATAACAAACTCCGGCC", "S000038"],
                    "D3GMAUX28": ["TATTTGCTTAA", "S000330"],
                    "CARG2ATAP3": ["CTTACCTTTCATGGATTA", "S000348"],
                    "OCETYPEIINTHISTONE": ["TCACGCGGATC", "S000268"],
                    "RSEPVGRP1": ["CAACTTTCATAT", "S000099"],
                    "GREGIONNTPRB1B": ["TGGCGGCTCTTATCTCACGTGATG", "S000364"],
                    "23BPZM27KDAZEIN": ["GACGTGTAAAGTAAATTTACAAC", "S000341"],
                    "GT1CONSENSUS": ["G[A|G][A|T]AA[A|T]", "S000198"],
                    "C2GMAUX28": ["AATAATAATAATAATAAATA", "S000327"],
                    "NRRBNEXTA": ["TAGTGGAT", "S000242"],
                    "SE1PVGRP18": ["ATAATGGGCCACACTGTGGGGCAT", "S000288"],
                    "MRNASTA1CRPSBD": ["CUCUUTGUTTUU", "S000274"],
                    "LRENPCABE": ["ACGTGGCA", "S000231"],
                    "D2GMAUX28": ["ATTTATATAAAT", "S000329"],
                    "ANAERO5CONSENSUS": ["TTCCCTGTT", "S000481"],
                    "MRNASTA2CRPSBD": ["UGAGUUG", "S000275"],
                    "TELOBOXATEEF1AA1": ["AAACCCTAA", "S000308"],
                    "L4DCPAL1": ["AATCTCCAACCA", "S000503"],
                    "ANAEROBICCISZMGAPC4": ["CGAAACCAGCAACGGTCCAG", "S000350"],
                    "PALBOXPPC": ["[C|T]T[C|T][C|T][A|C][A|C]C[A|C]A[A|C]C[A|C][A|C]C", "S000136"],
                    "EMHVCHORD": ["TGTAAAGT", "S000452"],
                    "ALF1NTPARC": ["TTACGCAAGCAATGACA", "S000238"],
                    "ANAERO1CONSENSUS": ["AAACAAA", "S000477"],
                    "D1GMAUX28": ["ACAGTTACTA", "S000328"],
                    "OSE2ROOTNODULE": ["CTCTT", "S000468"],
                    "OBF5ATGST6": ["ATCTTATGTCATTGATGACGACCTCC", "S000304"],
                    "CARG3ATAP3": ["CTTTCCATTTTTAGTAAC", "S000349"],
                    "LS7ATPR1": ["ACGTCATAGA", "S000322"],
                    "CONSERVED11NTZMATP1": ["ACGTATTAAAA", "S000301"],
                    "ABADESI1": ["[A|G]TACGTGGC[A|G]", "S000007"],
                    "ABADESI2": ["GGACGCGTGGC", "S000008"],
                    "LTRECOREATCOR15": ["CCGAC", "S000153"],
                    "-300ELEMENT": ["TG[A|C|T]AAA[A|G][G|T]", "S000122"],
                    "AMMORESVDCRNIA1": ["GGCCCCGGG", "S000376"],
                    "GT1CORE": ["GGTTAA", "S000125"],
                    "LEGUMINBOXLEGA5": ["TCCATAGCCATGCA[A|T][A|G]CTG[A|C]AGAATGTC", "S000060"],
                    "SUREAHVISO1": ["AAAACTAAGAAAGACCGATGGAAAA", "S000441"],
                    "ACGTSEED3": ["GTACGTGGCG", "S000019"],
                    "ACGTSEED2": ["ACACACGTCAA", "S000018"],
                    "NONAMERATH4": ["AGATCGACG", "S000147"],
                    "OCSENHANMOTIFAT": ["ACGTAAGCGCTTACGT", "S000074"],
                    "S1FBOXSORPS1L21": ["ATGGTA", "S000223"],
                    "CTRMCAMV35S": ["TCTCTCTCT", "S000460"],
                    "IRO2OS": ["CACGTGG", "S000505"],
                    "AGL2ATCONSENSUS": ["..[A|T].CCA[A|T][A|T][A|T][A|T]T[A|G]G[A|T][A|T]A.", "S000339"],
                    "SITE3SORPS1": ["AGTTAGTTAAAAGA", "S000212"],
                    "ELRE1PCPAL1": ["CTCCAACAAACCCCTTC", "S000306"],
                    "ABFOS": ["GCATCTTTACTTTAGCATC", "S000190"],
                    "CEREGLUBOX1PSLEGA": ["TGTTAAAGT", "S000032"],
                    "HDMOTIFPCPR2": ["CTAATTGTTTA", "S000233"],
                    "-300CORE": ["TGTAAAG", "S000001"],
                    "SV40COREENHAN": ["GTGG[A|T][A|T][A|C|T]G", "S000123"],
                    "ATHB5ATCORE": ["CAAT.ATTG", "S000371"],
                    "5659BOXLELAT5659": ["GAA[A|T]TTGTGA", "S000280"],
                    "ABRE2HVA22": ["CGCACGTGTC", "S000117"],
                    "GRAZMRAB28": ["CATGCCGCC", "S000220"],
                    "QARBNEXTA": ["AACGTGT", "S000244"],
                    "INRNTPSADB": ["[C|T]TCA.T[C|T][C|T]", "S000395"],
                    "BOXIINTPATPB": ["ATAGAA", "S000296"],
                    "PIIATGAPB": ["TTGGTTTTGATCAAAACCAA", "S000382"],
                    "BOX1PVCHS15": ["TAAAAGTTAAAAAC", "S000208"],
                    "SGBFGMGMAUX28": ["TCCACGTGTC", "S000287"],
                    "AGAMOUSATCONSENSUS": ["TT[A|G|T]CC[A|T][A|T][A|T][A|T][A|T][A|T]GG[A|C|T]AA", "S000342"],
                    "EECCRCAH1": ["GA.TT.C", "S000494"],
                    "TGACGTVMAMY": ["TGACGT", "S000377"],
                    "RAV1AAT": ["CAACA", "S000314"],
                    "SB3NPABC1": ["TTATGAACAGTAATT", "S000434"],
                    "BOX1PSGS2": ["ATAGAAATCAA", "S000222"],
                    "LEAFYATAG": ["CCAATGT", "S000432"],
                    "ABREMOTIFIIIOSRAB16B": ["GCCGCGTGGC", "S000291"],
                    "WBBOXPCWRKY1": ["TTTGAC[C|T]", "S000310"],
                    "ABREZMRAB28": ["CCACGTGG", "S000133"],
                    "UPRE1AT": ["ATTGGTCCACG", "S000428"],
                    "TATCCAYMOTIFOSRAMY3D": ["TATCCA[C|T]", "S000256"],
                    "DRE1COREZMRAB17": ["ACCGAGA", "S000401"],
                    "MARTBOX": ["TT[A|T]T[A|T]TT[A|T]TT", "S000067"],
                    "SP8BFIBSP8BIB": ["TACTATT", "S000184"],
                    "PROLAMINBOXOSGLUB1": ["TGCAAAG", "S000354"],
                    "ABRECE1HVA22": ["TGCCACCGG", "S000014"],
                    "ANAERO3CONSENSUS": ["TCATCAC", "S000479"],
                    "ABRERATCAL": ["[A|C]ACG[C|T]G[C|G|T]", "S000507"],
                    "SEF1MOTIF": ["ATATTTA[A|T][A|T]", "S000006"],
                    "AGL3ATCONSENSUS": ["TT[A|T]C[C|T]A[A|T][A|T][A|T][A|T]T[A|G]G[A|T]AA", "S000343"],
                    "RYREPEATBNNAPA": ["CATGCA", "S000264"],
                    "ERELEE4": ["A[A|T]TTCAAA", "S000037"],
                    "GBOXPC": ["ACCACGTGGC", "S000324"],
                    "OPAQUE2ZMB32": ["GATGA[C|T][A|G]TGG", "S000077"],
                    "UPRE2AT": ["CCACGTCATC", "S000429"],
                    "ABRECE3HVA1": ["ACGCGTGTCCTC", "S000141"],
                    "RHERPATEXPA7": ["[G|T]CACG[A|T]", "S000512"],
                    "MYBST1": ["GGATA", "S000180"],
                    "MYBPZM": ["CC[A|T]ACC", "S000179"],
                    "HSELIKENTGLN2": ["AGGAATTCCT", "S000057"],
                    "O2F3BE2S1": ["TCCACGTACT", "S000164"],
                    "IBOX": ["GATAAG", "S000124"],
                    "POLASIG3": ["AATAAT", "S000088"],
                    "POLASIG2": ["AATTAAA", "S000081"],
                    "POLASIG1": ["AATAAA", "S000080"],
                    "SE2PVGRP1": ["TT...GTAGCTAGTGTATTTGTAT", "S000101"],
                    "ACGTABREMOTIFAOSOSEM": ["TACGTGTC", "S000281"],
                    "HDZIPIIIAT": ["GTAAT[G|C]ATTAC", "S000475"],
                    "-10PEHVPSBD": ["TATTCT", "S000392"],
                    "MYBPLANT": ["[A|C]ACC[A|T]A[A|C]C", "S000167"],
                    "BOXINTPATPB": ["AATTCCATAGAATAGATAATA", "S000295"],
                    "CACGTGMOTIF": ["CACGTG", "S000042"],
                    "GARE4HVEPB1": ["GTAACAGAATGCTGG", "S000297"],
                    "CARGCW8GAT": ["C[A|T][A|T][A|T][A|T][A|T][A|T][A|T][A|T]G", "S000431"],
                    "SP8BFIBSP8AIB": ["ACTGTGTA", "S000183"],
                    "GBOXLERBCS": ["[A|C]CACGTGGC", "S000041"],
                    "HSRENTHSR203J": ["CAAAATTTTGTA", "S000466"],
                    "ABRETAEM": ["GGACACGTGGC", "S000015"],
                    "ANAERO2CONSENSUS": ["AGCAGC", "S000478"],
                    "AGCBOXNPGLB": ["AGCCGCC", "S000232"],
                    "SBOXATRBCS": ["CACCTCCA", "S000500"],
                    "RBCSBOX3PS": ["ATCATTTTCACT", "S000095"],
                    "INTRONUPPER": ["[A|C]AGGTAAGT", "S000085"],
                    "GLUTEBP1OS": ["AAGCAACACACAAC", "S000048"],
                    "ACGTROOT1": ["GCCACGTGGC", "S000016"],
                    "SEF4MOTIFGM7S": ["[A|G]TTTTT[A|G]", "S000103"],
                    "CDA1ATCAB2": ["CAAAACGC", "S000440"],
                    "AGL1ATCONSENSUS": [".TT[A|G|T]CC[A|T][A|T][A|T][A|T]..GG[A|T]AA.", "S000338"],
                    "E2FBNTRNR": ["GCGGCAAA", "S000367"],
                    "TATCCAOSAMY": ["TATCCA", "S000403"],
                    "WUSATAg": ["TTAATGG", "S000433"],
                    "LECPLEACS2": ["TAAAATAT", "S000465"],
                    "PREMOTIFNPCABE": ["ACCGGCCCACTT", "S000230"],
                    "NDEGMSAUR": ["CCATATGCCATGTCTCTCAATTGGTCCCAT", "S000359"],
                    "SARECAMV": ["CTGACGTAAGGGATGACGCAC", "S000156"],
                    "RBCSCONSENSUS": ["AATCCAA", "S000127"],
                    "COREOS": ["AA[G|T]AAT[A|T][C|T][A|G]TA[A|T]ATAAAA[A|C]TTTTAT[A|T]TA", "S000469"],
                    "REBETALGLHCB21": ["CGGATA", "S000363"],
                    "LREBOX3PSRBCS3": ["ACTATTTTCACTATC", "S000062"],
                    "IBOXCORENT": ["GATAAG[A|G]", "S000424"],
                    "PALINDROMICCBOXGM": ["TGACGTCA", "S000255"],
                    "DRE2COREZMRAB17": ["ACCGAC", "S000402"],
                    "SORLREP2AT": ["ATAAAACGT", "S000487"],
                    "JASE1ATOPR1": ["CGTCAATGAA", "S000388"],
                    "CATATGGMSAUR": ["CATATG", "S000370"],
                    "JASE2ATOPR1": ["CATACGTCGTCAA", "S000389"],
                    "ABASEED1": ["TGTTACGTGCC", "S000011"],
                    "GADOWNAT": ["ACGTGTC", "S000438"],
                    "MYBCOREATCYCB1": ["AACGG", "S000502"],
                    "PASNTPARA": ["TTACGCAAGCAATGACATCT", "S000336"],
                    "ABRE2HVA1": ["CCTACGTGGCGG", "S000134"],
                    "CACTFTPPCA1": ["[C|T]ACT", "S000449"],
                    "BOXCPSAS1": ["CTCCCAC", "S000226"],
                    "HBOXCONSENSUSPVCHS": ["CCTACC.......CT", "S000200"],
                    "PSREGIONZMZM13": ["TCGGCCACTATTTCTACGGGCAGCCAGACAAA", "S000253"],
                    "TBOXATGAPB": ["ACTTTG", "S000383"],
                    "MEJARELELOX": ["GATACA..AAT.TGATG", "S000151"],
                    "MARABOX1": ["AATAAA[C|T]AAA", "S000063"],
                    "ACEATCHS": ["GACACGTAGA", "S000355"],
                    "WBOXHVISO1": ["TGACT", "S000442"],
                    "ABREBZMRAB28": ["TCCACGTCTC", "S000219"],
                    "ABRE3OSRAB16": ["GTACGTGGCGC", "S000120"],
                    "3AF1BOXPSRBCS3": ["AAATAGATAAATAAAAACATT", "S000004"],
                    "GCN4OSGLUB1": ["TGAGTCA", "S000277"],
                    "WARBNEXTA": ["GTACGTGTTATAAAACGTGT", "S000241"],
                    "RNFG2OS": ["CCAGTGTGCCCCTGG", "S000189"],
                    "TATAPVTRNALEU": ["TTTATATA", "S000340"],
                    "BOX3PVCHS15": ["TATTGGTTACTAAA", "S000210"],
                    "ARE1": ["[A|G]GTGAC...GC", "S000022"],
                    "PALBOXLPC": ["[C|T]C[C|T][C|T]ACC[A|T]ACC", "S000138"],
                    "GARE2": ["[A|G]TAACA[A|G]A.TC[C|T]GG", "S000165"],
                    "HEXMOTIFTAH3H4": ["ACGTCA", "S000053"],
                    "WRKY71OS": ["TGAC", "S000447"],
                    "GLMHVCHORD": ["[A|G]TGA[G|C]TCAT", "S000451"],
                    "HSE": ["CT.GAA..TTC.AG", "S000054"],
                    "SORLREP3AT": ["TGTATATAT", "S000488"],
                    "REGION1OSOSEM": ["CGGCGGCCTCGCCACG", "S000300"],
                    "IDRSZMFER1": ["CACGAG[G|C]CC[G|T]CCAC", "S000445"],
                    "ANAERO4CONSENSUS": ["GTTT[A|C|T]GCAA", "S000480"],
                    "MYCCONSENSUSAT": ["CA..TG", "S000407"],
                    "23BPUASNSCYCB1": ["TTTATTTACCAAACGGTAACATC", "S000283"],
                    "EVENINGAT": ["AAAATATCT", "S000385"],
                    "GT1GMSCAM4": ["GAAAAA", "S000453"],
                    "AUXREPSIAA4": ["[G|T]GTCCCAT", "S000026"],
                    "PE1ASPHYA3": ["GAAATAGCAAATGTTAAAAATA", "S000196"],
                    "IBOXCORE": ["GATAA", "S000199"],
                    "PE3ASPHYA3": ["CAGCTCCCATGGCTCTCCCATCCGCGCCGGT", "S000197"],
                    "GARE1OSREP1": ["TAACAGA", "S000419"],
                    "ELRENTCHN50": ["GGTCA...AGTC", "S000155"],
                    "EIN3ATERF1": ["GGATTCAAGGGGCATGTATCTTGAATCC", "S000332"],
                    "ABREATRD22": ["[A|G][C|T]ACGTGG[C|T][A|G]", "S000013"],
                    "GMHDLGMVSPB": ["CATTAATTAG", "S000372"],
                    "AUXRETGA1GMGH3": ["TGACGTAA", "S000234"],
                    "CCAATBOX1": ["CCAAT", "S000030"],
                    "POLLEN2LELAT52": ["TCCACCATA", "S000246"],
                    "LTRE1HVBLT49": ["CCGAAA", "S000250"],
                    "GAGAGMGSA1": ["GAGAGAGAGAGAGAGAGA", "S000405"],
                    "-141NTG13": ["GCTTTTGATGACTTCAAACAC", "S000335"],
                    "SREATMSD": ["TTATCC", "S000470"],
                    "IDE2HVIDS2": ["TTGAACGGCAAGTTTCACGCTGTCACT", "S000464"],
                    "-314MOTIFZMSBE1": ["ACATAAAATAAAAAAAGGCA", "S000284"],
                    "-284MOTIFZMSBE1": ["CGTGCAAGCCCAAAGGCCAATCGGCCCAGA", "S000285"],
                    "P1BS": ["G.ATAT.C", "S000459"],
                    "CGTGTSPHZMC1": ["CGTGTCGTCCATGCAT", "S000294"],
                    "GBOXRELOSAMY3": ["CTACGTGGCCA", "S000206"],
                    "PRECONSCRHSP70A": ["[G|C]CGA[C|T].[A|G]...............[A|C|T][A|G|T]", "S000506"],
                    "ATHB2ATCONSENSUS": ["CAAT[G|C]ATTG", "S000318"],
                    "OCTAMERMOTIFTAH3H4": ["CGCGGATC", "S000076"],
                    "CGF1ATCAB2": ["GATAAAGATTACTTCAGATATAACAAACGTTAC", "S000213"],
                    "GT2OSPHYA": ["GCGGTAATT", "S000207"],
                    "BBOXSITE1STPAT": ["GCTAAACAAT", "S000398"],
                    "CCA1ATLHCB1": ["AA[A|C]AATCT", "S000149"],
                    "RNFG1OS": ["GATCATCGATC", "S000188"],
                    "HEXAMERATH4": ["CCGTCG", "S000146"],
                    "D4GMAUX28": ["TAGTGCTGT", "S000331"],
                    "ASF1MOTIFCAMV": ["TGACG", "S000024"],
                    "TATABOX3": ["TATTAAT", "S000110"],
                    "DREDR1ATRD29AB": ["TACCGACAT", "S000152"],
                    "MYB1AT": ["[A|T]AACCA", "S000408"],
                    "TATABOX4": ["TATATAA", "S000111"],
                    "HDZIP2ATATHB2": ["TAAT[A|C]ATTA", "S000373"],
                    "GATABOX": ["GATA", "S000039"],
                    "XYLAT": ["ACAAAGAA", "S000510"],
                    "ARELIKEGHPGDFR2": ["AGTTGAATGGGGGTGCA", "S000437"],
                    "NTBBF1ARROLB": ["ACTTTA", "S000273"],
                    "ABRE3HVA22": ["GCCACGTACA", "S000118"],
                    "PROLAMINBOX": ["CACATGTGTAAAGGT", "S000091"],
                    "ABREA2HVA1": ["CCTACGTGGC", "S000140"],
                    "O2F2BE2S1": ["GCCACCTCAT", "S000163"],
                    "NONAMERMOTIFTAH3H4": ["CATCCAACG", "S000071"],
                    "S2FSORPL21": ["CCATACATT", "S000166"],
                    "SORLIP1AT": ["GCCAC", "S000482"],
                    "ATHB6COREAT": ["CAATTATTA", "S000399"],
                    "SORLIP4AT": ["GTATGATGG", "S000485"],
                    "ABRECE3ZMRAB28": ["ACGCGCCTCCTC", "S000221"],
                    "E2FAT": ["T[C|T]TCCCGCC", "S000417"],
                    "IDE1HVIDS2": ["ATCAAGCATGCTTCTTGC", "S000463"],
                    "GLUTEBOX2OSGT3": ["CTTTTGTGTACCTTA", "S000129"],
                    "GLUTEBOX2OSGT2": ["TCCGTGTACCA", "S000047"],
                    "BOXLCOREDCPAL": ["ACC[A|T][A|T]CC", "S000492"],
                    "SB1NPABC1": ["CACTAACACAAAGTAA", "S000435"],
                    "-300MOTIFZMZEIN": ["[A|G]TGAGTCAT", "S000002"],
                    "WBOXNTERF3": ["TGAC[C|T]", "S000457"],
                    "MYCATERD1": ["CATGTG", "S000413"],
                    "TATABOXOSPAL": ["TATTTAA", "S000400"],
                    "IBOXLSCMCUCUMISIN": ["AGATATGATAAAA", "S000423"],
                    "ACGTCBOX": ["GACGTC", "S000131"],
                    "TRANSINITMONOCOTS": ["[A|G][A|C].AUGGC", "S000202"],
                    "CYTOSITECSHPRA": ["AAGATTGATTGAG", "S000261"],
                    "ASF1NTPARA": ["TTACGCAAGCAATGACAT", "S000240"],
                    "SEBFCONSSTPR10A": ["[C|T]TGTC[A|T]C", "S000391"],
                    "OBP1ATGST6": ["TACACTTTTGG", "S000305"],
                    "GLUTAACAOS": ["AACAAACTCTAT", "S000045"],
                    "DPBFCOREDCDC3": ["ACAC..G", "S000292"],
                    "TGTCACACMCUCUMISIN": ["TGTCACA", "S000422"],
                    "ATHB1ATCONSENSUS": ["CAAT[A|T]ATTG", "S000317"],
                    "RSRBNEXTA": ["CAAACTCGTATATCCAT", "S000243"],
                    "JERECRSTR": ["CTCTTAGACCGCCTTCTTTGAAAG", "S000384"],
                    "PROXBBNNAPA": ["CAAACACC", "S000263"],
                    "ROOTMOTIFTAPOX1": ["ATATT", "S000098"],
                    "SRENTTTO1": ["TGGTAGGTGAGAT", "S000271"],
                    "CAATBOX2": ["GGCCAATCT", "S000029"],
                    "CAATBOX1": ["CAAT", "S000028"],
                    "GCAACREPEATZMZEIN": ["GCAACGCAAC", "S000027"],
                    "UP1ATMSD": ["GGCCCA[A|T][A|T][A|T]", "S000471"],
                    "MYB2AT": ["TAACTG", "S000177"],
                    "SEF3MOTIFGM": ["AACCCA", "S000115"],
                    "OPAQUE2ZM22Z": ["TCCACGTAGA", "S000017"],
                    "PR2GCNT": ["TAA[A|G]AGCCGCC", "S000089"],
                    "NODCON2GM": ["CTCTT", "S000462"],
                    "TATCCACHVAL21": ["TATCCAC", "S000416"],
                    "LREBOXIPCCHS1": ["AACCTAACCT", "S000302"],
                    "AUXRETGA2GMGH3": ["TGACGTGGC", "S000235"],
                    "ZDNAFORMINGATCAB1": ["ATACGTGT", "S000321"],
                    "CEREGLUBOX2PSLEGA": ["TGAAAACT", "S000033"],
                    "ABREATCONSENSUS": ["[C|T]ACGTGGC", "S000406"],
                    "RE1ASPHYA3": ["CATGGGCGCGG", "S000195"],
                    "SITEIIBOSPCNA": ["TGGTCCCAC", "S000217"],
                    "RYREPEAT4": ["TCCATGCATGCAC", "S000010"],
                    "SITE1SORPS1": ["TCATGGTAACAATT", "S000211"],
                    "ACGTATERD1": ["ACGT", "S000415"],
                    "LS5ATPR1": ["TCTACGTCAC", "S000323"],
                    "RYREPEATLEGUMINBOX": ["CATGCA[C|T]", "S000100"],
                    "ELRECOREPCRP1": ["TTGACC", "S000142"],
                    "ACGTABOX": ["TACGTA", "S000130"],
                    "OCETYPEIIINTHISTONE": ["GATCCGCG..............ACCAATC[G|C]", "S000269"],
                    "DE1PSPRA2": ["GGATTTTACAGT", "S000379"],
                    "QELEMENTZMZM13": ["AGGTCA", "S000254"],
                    "TE2F2NTPCNA": ["ATTCCCGC", "S000397"],
                    "GARE2OSREP1": ["TAACGTA", "S000420"],
                    "LREBOX2PSRBCS3": ["TGTGTGGTTAATATG", "S000061"],
                    "BP5OSWX": ["CAACGTG", "S000436"],
                    "UP2ATMSD": ["AAACCCTA", "S000472"],
                    "B2GMAUX28": ["CTTGTCGTCA", "S000325"],
                    "CARGNCAT": ["CC[A|T][A|T][A|T][A|T][A|T][A|T][A|T][A|T]GG", "S000446"],
                    "EMBP1TAEM": ["CACGTGGC", "S000119"],
                    "GT1MOTIFPSRBCS": ["[G|T][A|T]GTG[A|G][A|T]AA[A|T][A|G][A|T]", "S000051"],
                    "27BPDRCONSENSUSPS25S": ["TCCGCC[A|T]CTTGTATTCGTTGCGTTG[A|C]A", "S000286"],
                    "CPRFPCCHS": ["CCACGTGGCC", "S000313"],
                    "C1GMAUX28": ["TGAAAACAGTGAGTTA", "S000326"],
                    "TATABOX5": ["TTATTT", "S000203"],
                    "AGTACSAO": ["AAAAAGTAAAAAGTAAAAAAGTAAAAAG", "S000258"],
                    "AACAOSGLUB1": ["CAACAAACTATATC", "S000276"],
                    "GLUTECOREOS": ["CTTTCGTGTAC", "S000050"],
                    "SURE2STPAT21": ["AATACTAAT", "S000185"],
                    "LBOXLERBCS": ["AAATTAACCAA", "S000126"],
                    "T/GBOXATPIN2": ["AACGTG", "S000458"],
                    "TDBA12NTCHN50": ["TGACTTTCTGAC", "S000266"],
                    "EBOXBNNAPA": ["CA..TG", "S000144"],
                    "MYBCORE": ["C.GTT[A|G]", "S000176"],
                    "OCSELEMENTAT": ["TGACG[C|T]AAG[G|C][A|G][A|C]T[G|T]ACG[C|T][A|C][A|C]", "S000158"],
                    "DR5GMGH3": ["CCTTTTGTCTC", "S000337"],
                    "TL1ATSAR": ["CTGAAGAAGAA", "S000473"],
                    "ABAREG2": ["ATGTACGAAGC", "S000009"],
                    "TCA1MOTIF": ["TCATCTTCTT", "S000159"],
                    "AS1CAMV": ["CCACTGACGTAAGGGATGACGCACAATCC", "S000023"],
                    "INTRONLOWER": ["TGCAGG", "S000086"],
                    "CIACADIANLELHC": ["CAA....ATC", "S000252"],
                    "GAGA8HVBKN3": ["GAGAGAGAGAGAGAGA", "S000427"],
                    "SORLIP2AT": ["GGGCC", "S000483"],
                    "BOXC'PSAS1": ["TCCCGGTACACACTTCTT", "S000227"],
                    "RBCSGBOXPS": ["CACATGGCACT", "S000096"],
                    "PALBOXAPC": ["CCGTCC", "S000137"],
                    "CARGATCONSENSUS": ["CC[A|T][A|T][A|T][A|T][A|T][A|T]GG", "S000404"],
                    "RYREPEATGMGY2": ["CATGCAT", "S000105"],
                    "ABREAZMRAB28": ["GCCACGTGGG", "S000218"],
                    "ABREOSRAB21": ["ACGT[G|C][G|C][G|C]C", "S000012"],
                    "TRANSINITDICOTS": ["A[A|C].AUGGC", "S000201"],
                    "ARR1AT": [".GATT", "S000454"],
                    "O2F1BE2S1": ["TCCACGTCGA", "S000162"],
                    "ARFAT": ["TGTCTC", "S000270"],
                    "SITEIIAOSPCNA": ["TGGGCCCGT", "S000216"],
                    "ARECOREZMGAPC4": ["AGCAACGGTC", "S000393"],
                    "LREBOXIIPCCHS1": ["TCCACGTGGC", "S000303"],
                    "MYBGAHV": ["TAACAAA", "S000181"],
                    "WRECSAA01": ["AA[A|T]GTATC[G|C]A", "S000496"],
                    "E2FCONSENSUS": ["[A|T]TT[G|C][G|C]C[G|C][G|C]", "S000476"],
                    "BS1EGCCR": ["AGCGGG", "S000352"],
                    "DRECRTCOREAT": ["[A|G]CCGAC", "S000418"],
                    "YREGIONNTPRB1B": ["TGTGACATTGAAATTCTTTGACTTTA", "S000365"],
                    "MYB2CONSENSUSAT": ["[C|T]AAC[G|T]G", "S000409"],
                    "MNF1ZMPPC1": ["GTGCCCTT", "S000251"],
                    "CE3OSOSEM": ["AACGCGTGTC", "S000282"],
                    "MYBATRD22": ["CTAACCA", "S000175"],
                    "OCSGMGH24": ["CGGTTTACGTAATCTCTTACATCA", "S000346"],
                    "PREATPRODH": ["ACTCAT", "S000450"],
                    "ABRE3HVA1": ["GCAACGTGTC", "S000135"],
                    "ACIIIPVPAL2": ["GTTAGGTTC", "S000194"],
                    "WBOXGACAD1A": ["AGTCAAAATTGACC", "S000448"],
                    "AMMORESIVDCRNIA1": ["CGAACTT", "S000375"],
                    "L1BOXATPDF1": ["TAAATG[C|T]A", "S000386"],
                    "ACGTOSGLUB1": ["GTACGTG", "S000278"],
                    "ACGTTBOX": ["AACGTT", "S000132"],
                    "ALF2NTPARB": ["TGAGGAGACTTGTGAGGT", "S000239"],
                    "ASF1ATNOS": ["TGAGCTAAGCACATACGTCAG", "S000073"],
                    "ELEMENT2GMLBC3": ["CTTAAATTATTTATTT", "S000320"],
                    "20NTNTNOS": ["TGAGCTAAGCACATACGTCA", "S000312"],
                    "RBCSBOX2PS": ["GTGTGGTTAATATG", "S000094"],
                    "BOXBPSAS1": ["AAACGACACCGTTT", "S000225"],
                    "ABREMOTIFAOSOSEM": ["TACGTGTC", "S000299"],
                    "GBOXSORBCS1": ["TCCACGTGGT", "S000380"],
                    "CGACGOSAMY3": ["CGACG", "S000205"],
                    "BOXICHS": ["GTCC[A|C]TC[A|C]AACCTA[A|C]C", "S000228"],
                    "GLUTEBP2OS": ["ATGCTCAATAGATATAAGT", "S000049"],
                    "TEF1BOXATA1": ["ACAGGGGCATAATGGTAATTTAAA", "S000311"],
                    "VSF1PVGRP18": ["GCTCCGTTG", "S000249"],
                    "PYRIMIDINEBOXOSRAMY1A": ["CCTTTT", "S000259"],
                    "AMYBOX1": ["TAACA[A|G]A", "S000020"],
                    "AMYBOX2": ["TATCCAT", "S000021"],
                    "ACIPVPAL2": ["CCCACCTACC", "S000192"],
                    "GGTCCCATGMSAUR": ["GGTCCCAT", "S000360"],
                    "SPHCOREZMC1": ["TCCATGCAT", "S000154"],
                    "AS1LIKECSHPRA": ["AAATGACGAAAATGC", "S000260"],
                    "LTREATLTI78": ["ACCGACA", "S000157"],
                    "TOPOISOM": ["GT.[A|T]A[C|T]ATT.AT..G", "S000112"],
                    "CGCGBOXAT": ["[A|C|G]CGCG[C|G|T]", "S000501"],
                    "ESPASGL01": ["ACATGTCATCATGT", "S000509"],
                    "GBOX10NT": ["GCCACGTGCC", "S000272"],
                    "UPRMOTIFIAT": ["CCACGTCA", "S000425"],
                    "CBFHV": ["[A|G][C|T]CGAC", "S000497"],
                    "ABREBNNAPA": ["CGCCACGTGTCC", "S000145"],
                    "RAV1BAT": ["CACCTG", "S000315"],
                    "ELRE2PCPAL1": ["ATTCTCACCTACCA", "S000307"],
                    "CPBCSPOR": ["TATTAG", "S000491"],
                    "OSE1ROOTNODULE": ["AAAGAT", "S000467"],
                    "NAPINMOTIFBN": ["TACACAT", "S000070"],
                    "OCTAMOTIF2": ["CGCGGCAT", "S000116"],
                    "OCETYPEINTHISTONE": ["CCACGTCA.CGATCCGCG", "S000267"],
                    "CANBNNAPA": ["C.AACAC", "S000148"],
                    "ABREDISTBBNNAPA": ["GCCACTTGTC", "S000262"],
                    "2SSEEDPROTBANAPA": ["CAAACAC", "S000143"],
                    "REALPHALGLHCB21": ["AACCAA", "S000362"],
                    "CRTDREHVCBF2": ["GTCGAC", "S000411"],
                    "TEFBOXATEEF1AA1": ["AGGGGCATAATGGTAA", "S000309"],
                    "AAGACGTAGATACL12": ["AAGACGTAG", "S000344"],
                    "SITEIOSPCNA": ["CCAGGTGG", "S000224"],
                    "BOXIIPCCHS": ["ACGTGGC", "S000229"],
                    "TGA1ANTPR1A": ["CGTCATCGAGATGACG", "S000247"],
                    "MYB26PS": ["GTTAGGTT", "S000182"]}

    with open(seeder_result_file, 'r') as seeder_result:
        PLACEd_seeder_filename = str(seeder_result_file) + ".html"
        PLACEd_seeder_file = open(PLACEd_seeder_filename, 'w')

        head = """
    <!DOCTYPE html>
        <html>
    <head>
        <title>Seeder results</title>
    </head>

    <body>
        <h1>Your Seeder results:</h1>
        <DIV style="font-family: 'Courier', 'monospace'; font-size:14px">\n
    """
        tail = """
    </DIV>
        </body>
        </html>
        """

        PLACEd_seeder_file.write(head)

        for line in seeder_result:
            line = line.strip().split()  # motif would be line[1]
            motif = line[1]
            motif_revcomp = str(Seq(motif).reverse_complement())
            line_to_write = '\t'.join(line)

            for m in PLACE_db:
                place_link = ''
                place_revcomp_link = ''
                for result in re.finditer(PLACE_db[m][0], motif):
                    place_link = str('<a href="http://www.dna.affrc.go.jp/sigscan/disp.cgi?%s">%s</a> ' % (PLACE_db[m][1], m)) + '\t' + result.group(0)
                for result in re.finditer(PLACE_db[m][0], motif_revcomp):
                    place_revcomp_link = str('<a href="http://www.dna.affrc.go.jp/sigscan/disp.cgi?%s">%s</a>' % (PLACE_db[m][1], m)) + '\t' + result.group(0) + '\t' + "Note: reverse complement"

                if place_link != '' or place_revcomp_link != '':
                    line_to_write_final = line_to_write + '\t' + place_link + \
                        '\t' + place_revcomp_link + '<br />' + '\n'

                    PLACEd_seeder_file.write(line_to_write_final)

        PLACEd_seeder_file.write(tail)
        PLACEd_seeder_file.close()

seeder_result_analysis("TESTFILE1.txt")
