MIN_JC = 0.5
FREQ_WINDOW = 0.10
CD_HIT_PATH = '/home/marcelo/Documentos/Programas/CDHIT/cdhit-master/./cd-hit'
MIN_CORR_FREQ = 0.05
MIN_CORR_COUNT = 100

properties = {}
properties[0] = "A"
properties[1] = "C"
properties[2] = "D"
properties[3] = "E"
properties[4] = "F"
properties[5] = "G"
properties[6] = "H"
properties[7] = "I"
properties[8] = "K"
properties[9] = "L"
properties[10] = "M"
properties[11] = "N"
properties[12] = "P"
properties[13] = "Q"
properties[14] = "R"
properties[15] = "S"
properties[16] = "T"
properties[17] = "V"
properties[18] = "W"
properties[19] = "Y"
properties[20] = "Amide"
properties[21] = "Aliphatic"
properties[22] = "Basic"
properties[23] = "Hydroxyl"
properties[24] = "Sulfur"
properties[25] = "Non-Polar"
properties[26] = "Polar"
properties[27] = "Hydrophobic"
properties[28] = "Hydrophilic"
properties[29] = "Pos.Charged"
properties[30] = "Neg.Charged"
properties[31] = "VerySmall"
properties[32] = "Small"
properties[33] = "MediumA"
properties[34] = "MediumB"
properties[35] = "Aromatic"
properties[36] = "ND"
properties[37] = "QE"

setsIdList = {}
setsIdList["A"] = [0,21,25,31]
setsIdList["C"] = [1,24,26,32]
setsIdList["D"] = [2,28,30,32,36]
setsIdList["E"] = [3,30,33,37]
setsIdList["F"] = [4,25,27,35]
setsIdList["G"] = [5,21,25,31]
setsIdList["H"] = [6,22,33]
setsIdList["I"] = [7,25,27,34]
setsIdList["K"] = [8,22,28,29,34]
setsIdList["L"] = [9,21,25,27,34]
setsIdList["M"] = [10,24,25,27,34]
setsIdList["N"] = [11,20,26,28,32,36]
setsIdList["P"] = [12,25,28,32]
setsIdList["Q"] = [13,20,26,28,33,37]
setsIdList["R"] = [14,22,28,29,34]
setsIdList["S"] = [15,23,26,31]
setsIdList["T"] = [16,23,26,32]
setsIdList["V"] = [17,21,25,27,33]
setsIdList["W"] = [18,25,27,35]
setsIdList["Y"] = [19,21,23,26,35]

setsAA = {}
setsAA[0] = "A"
setsAA[1] = "C"
setsAA[2] = "D"
setsAA[3] = "E"
setsAA[4] = "F"
setsAA[5] = "G"
setsAA[6] = "H"
setsAA[7] = "I"
setsAA[8] = "K"
setsAA[9] = "L"
setsAA[10] = "M"
setsAA[11] = "N"
setsAA[12] = "P"
setsAA[13] = "Q"
setsAA[14] = "R"
setsAA[15] = "S"
setsAA[16] = "T"
setsAA[17] = "V"
setsAA[18] = "W"
setsAA[19] = "Y"
setsAA[20] = "NQ"
setsAA[21] = "GAVLY"
setsAA[22] = "HKR"
setsAA[23] = "STY"
setsAA[24] = "CM"
setsAA[25] = "FGVLAIPMW"
setsAA[26] = "YSNTQC"
setsAA[27] = "LIFWVM"
setsAA[28] = "RKNQPD"
setsAA[29] = "KR"
setsAA[30] = "DE"
setsAA[31] = "GAS"
setsAA[32] = "CDNPT"
setsAA[33] = "EVQH"
setsAA[34] = "MILKR"
setsAA[35] = "FYW"
setsAA[36] = "ND"
setsAA[37] = "QE"

alphabet = {}
alphabet["A"] = 1
alphabet["C"] = 1
alphabet["D"] = 1
alphabet["E"] = 1
alphabet["F"] = 1
alphabet["G"] = 1
alphabet["H"] = 1
alphabet["I"] = 1
alphabet["K"] = 1
alphabet["L"] = 1
alphabet["M"] = 1
alphabet["N"] = 1
alphabet["P"] = 1
alphabet["Q"] = 1
alphabet["R"] = 1
alphabet["S"] = 1
alphabet["T"] = 1
alphabet["V"] = 1
alphabet["W"] = 1
alphabet["Y"] = 1
alphabet["Amide"] = 2
alphabet["Aliphatic"] = 5
alphabet["Basic"] = 3
alphabet["Hydroxyl"] = 3
alphabet["Sulfur"] = 2
alphabet["Non-Polar"] = 9
alphabet["Polar"] = 6.1
alphabet["Hydrophobic"] = 6
alphabet["Hydrophilic"] = 6
alphabet["Pos.Charged"] = 2
alphabet["Neg.Charged"] = 2
alphabet["VerySmall"] = 3.2
alphabet["Small"] = 5.1
alphabet["MediumA"] = 4
alphabet["MediumB"] = 5.1
alphabet["Aromatic"] = 3.1
alphabet["ND"] = 2.2
alphabet["QE"] = 2.2