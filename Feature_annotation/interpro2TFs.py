#!/usr/bin/python

'''
This program searches introproscan results tables and extyracts those genes with
functional annotations associated with transcription factor domains.
'''


#-----------------------------------------------------
# Step 1
# Import variables & load input files
#-----------------------------------------------------

import sys
import argparse
import re
from sets import Set
from collections import defaultdict
from operator import itemgetter

ap = argparse.ArgumentParser()

ap.add_argument('--InterPro',required=True,type=str,help='The Interproscan functional annotation .tsv file')

conf = ap.parse_args()


with open(conf.InterPro) as f:
    InterPro_lines = f.readlines()

column_list=[]


#-----------------------------------------------------
# Step 2
#
#-----------------------------------------------------
# This uses superfamily and PFAM domains from DBD and
# IPR terms from Shelest 2017
# These terms were noted to miss the IPR terms for
# Fungal-specific transcription factors IPR007219
# and IPR021858, which have been included.

TF_list = [
('IPR007219', 'Fungal-specific transcription factor'),
('IPR021858', 'Fungal-specific transcription factor'),
('SSF101941', 'NAC domain'),
('SSF103612', 'SBT domain'),
('SSF103637', 'CCHHC domain'),
('SSF109709', 'KorB DNA-binding domain-like'),
('SSF46689', 'Homeodomain-like'),
('SSF46785', '"Winged helix" DNA-binding domain'),
('SSF46894', 'C-terminal effector domain of the bipartite response regulators'),
('SSF46955', 'Putative DNA-binding domain'),
('SSF47413', 'lambda repressor-like DNA-binding domains'),
('SSF47454', 'A DNA-binding domain in eukaryotic transcription factors'),
('SSF47459', 'HLH, helix-loop-helix DNA-binding domain'),
('SSF47598', 'Ribbon-helix-helix'),
('SSF47655', 'STAT'),
('SSF48295', 'TrpR-like'),
('SSF49417', 'p53-like transcription factors'),
('SSF50249', 'Nucleic acid-binding proteins'),
('SSF51215', 'Regulatory protein AraC'),
('SSF54171', 'DNA-binding domain'),
('SSF54447', 'ssDNA-binding transcriptional regulator domain'),
('SSF54518', 'Transcriptional factor tubby, C-terminal domain'),
('SSF54616', 'DNA-binding domain of Mlu1-box binding protein MBP1'),
('SSF54957', 'Viral DNA-binding domain'),
('SSF55021', 'ACT-like'),
('SSF55455', 'SRF-like'),
('SSF56366', 'SMAD MH1 domain'),
('SSF56548', 'Conserved core of transcriptional regulatory protein vp16'),
('SSF57667', 'C2H2 and C2HC zinc fingers'),
('SSF57701', 'Zn2/Cys6 DNA-binding domain'),
('SSF57716', 'Glucocorticoid receptor-like (DNA-binding domain)'),
('SSF57879', 'Zinc domain conserved in yeast copper-regulated transcription factors'),
('SSF63592', 'Flagellar transcriptional activator FlhD'),
('SSF63763', 'SAND domain-like'),
('SSF68989', 'Hemolysin expression modulating protein HHA'),
('SSF69652', 'DNA-binding C-terminal domain of the transcription factor MotA'),
('SSF82927', 'Cysteine-rich DNA binding domain, (DM domain)'),
('SSF89447', 'AbrB/MazE/MraZ-like'),
('SSF89915', 'DNA-binding protein Tfx'),
('SSF90073', 'GCM domain'),
('PF00010', 'Helix-loop-helix DNA-binding domain'),
('PF00046', 'Homeobox domain'),
('PF00096', 'Zinc finger, C2H2 type'),
('PF00105', 'Zinc finger, C4 type (two domains)'),
('PF00126', 'Bacterial regulatory helix-turn-helix protein, lysR family'),
('PF00157', 'Pou domain - N-terminal to homeobox domain'),
('PF00165', 'Bacterial regulatory helix-turn-helix proteins, AraC family'),
('PF00170', 'bZIP transcription factor'),
('PF00178', 'Ets-domain'),
('PF00196', 'Bacterial regulatory proteins, luxR family'),
('PF00250', 'Fork head domain'),
('PF00292', '\'Paired box\' domain'),
('PF00313', '\'Cold-shock\' DNA-binding domain'),
('PF00319', 'SRF-type transcription factor (DNA-binding and dimerisation domain)'),
('PF00320', 'GATA zinc finger'),
('PF00325', 'Bacterial regulatory proteins, crp family'),
('PF00356', 'Bacterial regulatory proteins, lacI family'),
('PF00376', 'MerR family regulatory protein'),
('PF00392', 'Bacterial regulatory proteins, gntR family'),
('PF00447', 'HSF-type DNA-binding'),
('PF00486', 'Transcriptional regulatory protein, C terminal'),
('PF00554', 'Rel homology domain (RHD)'),
('PF00605', 'Interferon regulatory factor transcription factor'),
('PF00649', 'Copper fist DNA binding domain'),
('PF00751', 'DM DNA binding domain'),
('PF00847', 'AP2 domain'),
('PF00853', 'Runt domain'),
('PF00870', 'P53 DNA-binding domain'),
('PF00907', 'T-box'),
('PF01022', 'Bacterial regulatory protein, arsR family'),
('PF01047', 'MarR family'),
('PF01056', 'Myc amino-terminal region'),
('PF01167', 'Tub family'),
('PF01285', 'TEA/ATTS domain family'),
('PF01316', 'Arginine repressor, DNA binding domain'),
('PF01325', 'Iron dependent repressor, N-terminal DNA binding domain'),
('PF01340', 'Met Apo-repressor, MetJ'),
('PF01342', 'SAND domain'),
('PF01371', 'Trp repressor protein'),
('PF01381', 'Helix-turn-helix'),
('PF01418', 'Helix-turn-helix domain, rpiR family'),
('PF01422', 'NF-X1 type zinc finger'),
('PF01475', 'Ferric uptake regulator family'),
('PF01530', 'Zinc finger, C2HC type'),
('PF01586', 'Myogenic Basic domain'),
('PF01638', 'HxlR-like helix-turn-helix'),
('PF01726', 'LexA DNA binding domain'),
('PF01754', 'A20-like zinc finger'),
('PF01978', 'Sugar-specific transcriptional regulator TrmB'),
('PF02045', 'CCAAT-binding transcription factor (CBF-B/NF-YA) subunit B'),
('PF02082', 'Transcriptional regulator'),
('PF02135', 'TAZ zinc finger'),
('PF02232', 'Alpha trans-inducing protein (Alpha-TIF)'),
('PF02236', 'Viral DNA-binding protein, all alpha domain'),
('PF02257', 'RFX DNA-binding domain'),
('PF02292', 'APSES domain'),
('PF02319', 'E2F/DP family winged-helix DNA-binding domain'),
('PF02362', 'B3 DNA binding domain'),
('PF02365', 'No apical meristem (NAM) protein'),
('PF02376', 'CUT domain'),
('PF02456', 'Adenovirus IVa2 protein'),
('PF02467', 'Transcription factor WhiB'),
('PF02701', 'Dof domain, zinc finger'),
('PF02703', 'Early E1A protein'),
('PF02791', 'DDT domain'),
('PF02864', 'STAT protein, DNA binding domain'),
('PF02891', 'MIZ zinc finger'),
('PF02892', 'BED zinc finger'),
('PF02928', 'C5HC2 zinc finger'),
('PF02944', 'BESS motif'),
('PF02954', 'Bacterial regulatory protein, Fis family'),
('PF03041', 'lef-2'),
('PF03106', 'WRKY DNA -binding domain'),
('PF03110', 'SBP domain'),
('PF03131', 'bZIP Maf transcription factor'),
('PF03299', 'Transcription factor AP-2'),
('PF03333', 'Adhesin biosynthesis transcription regulatory protein'),
('PF03343', 'SART-1 family'),
('PF03551', 'Transcriptional regulator PadR-like family'),
('PF03584', 'Herpesvirus ICP4-like protein N-terminal region'),
('PF03615', 'GCM motif protein'),
('PF03634', 'TCP family transcription factor'),
('PF03749', 'Sugar fermentation stimulation protein'),
('PF03859', 'CG-1 domain'),
('PF03869', 'Arc-like DNA binding domain'),
('PF03965', 'Penicillinase repressor'),
('PF04014', 'SpoVT / AbrB like domain'),
('PF04024', 'PspC domain'),
('PF04054', 'CCR4-Not complex component, Not1'),
('PF04247', 'Invasion gene expression up-regulator, SirB'),
('PF04299', 'Putative FMN-binding domain'),
('PF04353', 'Regulator of RNA polymerase sigma(70) subunit, Rsd/AlgQ'),
('PF04383', 'KilA-N domain'),
('PF04397', 'LytTr DNA-binding domain'),
('PF04441', 'Poxvirus early transcription factor (VETF), large subunit'),
('PF04516', 'CP2 transcription factor'),
('PF04606', 'Phage transcriptional activator, Ogr/Delta'),
('PF04689', 'DNA binding protein S1FA'),
('PF04690', 'YABBY protein'),
('PF04745', 'VITF-3 subunit protein'),
('PF04761', 'Lactococcus bacteriophage putative transcription regulator'),
('PF04769', 'Mating-type protein MAT alpha 1'),
('PF04838', 'Baculoviridae late expression factor 5'),
('PF04873', 'Ethylene insensitive 3'),
('PF04947', 'Poxvirus Late Transcription Factor VLTF3 like'),
('PF04967', 'HTH DNA binding domain'),
('PF05043', 'Mga helix-turn-helix domain'),
('PF05044', 'Homeobox prospero-like protein (PROX1)'),
('PF05110', 'AF-4 proto-oncoprotein'),
('PF05224', 'NDT80 / PhoG like DNA-binding  family'),
('PF05225', 'helix-turn-helix, Psq domain'),
('PF05234', 'UAF complex subunit Rrn10'),
('PF05247', 'Flagellar transcriptional activator (FlhD)'),
('PF05280', 'Flagellar transcriptional activator (FlhC)'),
('PF05290', 'Baculovirus immediate-early protein (IE-0)'),
('PF05311', 'Baculovirus 33KDa late protein (PP31)'),
('PF05321', 'Haemolysin expression modulating protein'),
('PF05443', 'ROS/MUCR transcriptional regulator protein'),
('PF05459', 'Herpesvirus transcriptional regulator family'),
('PF05464', 'Phi-29-like late genes activator (early protein GP4)'),
('PF05718', 'Poxvirus intermediate transcription factor'),
('PF05764', 'YL1 nuclear protein'),
('PF05848', 'Firmicute transcriptional repressor of class III stress genes (CtsR)'),
('PF06018', 'CodY GAF-like domain'),
('PF06054', 'Competence protein CoiA-like family'),
('PF06069', 'PerC transcriptional activator'),
('PF06116', 'Transcriptional activator RinB'),
('PF06200', 'ZIM motif'),
('PF06320', 'GCN5-like protein 1 (GCN5L1)'),
('PF06338', 'ComK protein'),
('PF06573', 'Churchill protein'),
('PF06719', 'AraC-type transcriptional regulator N-terminus'),
('PF06818', 'Fez1'),
('PF06839', 'GRF zinc finger'),
('PF06943', 'LSD1 zinc finger'),
('PF07093', 'SGT1 protein'),
('PF07417', 'Transcriptional regulator Crl'),
('PF07704', 'Rv0623-like transcription factor'),
('PF07716', 'Basic region leucine zipper'),
('PF07750', 'GcrA cell cycle regulator'),
('PF07764', 'Omega Transcriptional Repressor'),
('PF07879', 'PHB/PHA accumulation regulator DNA-binding domain'),
('PF08220', 'DeoR-like helix-turn-helix domain'),
('PF08222', 'CodY helix-turn-helix domain'),
('PF08270', 'M protein trans-acting positive regulator (MGA) PRD domain'),
('PF08279', 'HTH domain'),
('PF08280', 'M protein trans-acting positive regulator (MGA) HTH domain'),
('IPR000005', 'Helix-turn-helix, AraC type'),
('IPR000007', 'Tubby, C-terminal'),
('IPR000197', 'Zinc finger, TAZ-type'),
('IPR000232', 'Heat shock factor (HSF)-type, DNA-binding'),
('IPR000327', 'POU-specific'),
('IPR000418', 'Ets'),
('IPR000551', 'Bacterial regulatory protein, MerR'),
('IPR000571', 'Zinc finger CCCH-type'),
('IPR000679', 'Zinc finger, GATA-type'),
('IPR000792', 'Bacterial regulatory protein, LuxR'),
('IPR000814', 'TATA-box binding'),
('IPR000818', 'TEA/ATTS'),
('IPR000835', 'Bacterial regulatory protein, MarR'),
('IPR000843', 'Bacterial regulatory protein, LacI'),
('IPR000944', 'Transcriptional regulator, Rrf2'),
('IPR000967', 'Zinc finger, NF-X1-type'),
('IPR001034', 'Bacterial regulatory protein, DeoR N-terminal'),
('IPR001083', 'Copper fist DNA-binding*'),
('IPR001138', 'Zn2 Cys6 Zn_cluster*'),
('IPR001275', 'DM DNA-binding'),
('IPR001289', 'CCAAT-binding TF, subunit B'),
('IPR001356', 'Homeobox'),
('IPR001387', 'Helix-turn-helix type 3'),
('IPR001471', 'Pathogenesis-related TF and ERF, DBD'),
('IPR001523', 'Paired box protein, N-terminal'),
('IPR001699', 'Transcription factor, T-box'),
('IPR001766', 'Fork head transcription factor'),
('IPR001808', 'Bacterial regulatory protein, Crp'),
('IPR001845', 'Bacterial regulatory protein, ArsR'),
('IPR001878', 'Zinc finger, CCHC-type'),
('IPR002059', 'Cold-shock protein, DNA-binding'),
('IPR002100', 'Transcription factor, MADS-box'),
('IPR002197', 'Helix-turn-helix, Fis-type'),
('IPR002653', 'Zinc finger, A20-type'),
('IPR003150', 'DNA-binding RFX'),
('IPR003163', 'APSES-type DNA-binding domain*'),
('IPR003316', 'E2F/dimerisation partner (TDP)'),
('IPR003656', 'Zinc finger, BED-type predicted'),
('IPR003657', 'DNA-binding WRKY'),
('IPR003902', 'Transcriptional regulator, GCM-like'),
('IPR003958', 'TF CBF/NF-Y/archaeal histone'),
('IPR004022', 'DDT'),
('IPR004181', 'Zinc finger, MIZ-type'),
('IPR004198', 'Zinc finger, C5HC2-type'),
('IPR004333', 'Transcription factor, SBP-box'),
('IPR004645', 'DNA-binding protein Tfx'),
('IPR004823', 'TATA box binding protein associated factor (TAF)'),
('IPR004826', 'Maf transcription factor'),
('IPR004827', 'Basic-leucine zipper (bZIP) TF'),
('IPR005011', 'SART-1 protein'),
('IPR006780', 'YABBY protein'),
('IPR006856', 'Mating-type protein MAT alpha 1*'),
('IPR007087', 'Zinc finger, C2H2-type'),
('IPR007196', 'CCR4-Not complex component, Not1'),
('IPR007396', 'Negative transcriptional regulator'),
('IPR007604', 'CP2 transcription factor'),
('IPR007889', 'Helix-turn-helix, Psq'),
('IPR008895', 'YL1 nuclear'),
('IPR008917', 'Eukaryotic transcription factor, Skn-1-like'),
('IPR008967', 'p53-like transcription factor, DNA-binding'),
('IPR009044', 'ssDNA-binding transcriptional regulator'),
('IPR009057', 'Homeodomain-like'),
('IPR009061', 'Putative DNA binding'),
('IPR009395', 'GCN5-like 1'),
('IPR010666', 'Zinc finger, GRF-type'),
('IPR010770', 'SGT1'),
('IPR010919', 'SAND-like'),
('IPR010921', 'Trp repressor/replication initiator'),
('IPR010982', 'Lambda repressor-like, DNA-binding'),
('IPR010985', 'Ribbon-helix-helix'),
('IPR011598', 'Helix-loop-helix DNA-binding'),
('IPR012294', 'Transcription factor TFIID, C-terminal'),
('IPR013921', 'TATA-binding related factor'),
('IPR013932', 'TATA-binding protein interacting (TIP20)'),
('IPR015988', 'STAT transcription factor, coiled coil'),
('IPR016032', 'Signal transduction response regulator, C-term. effector'),
('IPR016177', 'DNA-binding, integrase-type'),
('IPR024061', 'NDT80 DNA-binding domain'),
('IPR025659', 'Tubby C-terminal-like domain')
]

TF_dict = defaultdict(list)
for k, entry in TF_list:
    TF_dict[k].append(entry)

#-----------------------------------------------------
# Step 3
#
#-----------------------------------------------------

outline_set = set()
outlines = []
for line in InterPro_lines:
    line = line.rstrip()
    split_line = line.split()
    gene_id = split_line[0]
    for col in split_line[4:]:
        if 'IPR' in col or 'PFAM' in col or 'SSF' in col:
            if TF_dict[col]:
                pred_func = "".join(TF_dict[col])
                outline = "\t".join([gene_id, col, pred_func])
                if outline in outline_set:
                    continue
                else:
                    outline_set.add(outline)
                    outlines.append(outline)

# outlines = sorted(set(outlines))
print "\n".join(outlines)
