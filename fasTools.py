#! /usr/bin/python
# -*-coding:Utf-8 -*

# ---------------------------------------------------------------------------- #
#
# Packages regroupant différentes fonctions pour manipuler des séquences
# nucléotidiques au format fasta
#
# ---------------------------------------------------------------------------- #
#
# Auteur: Thomas Goeury
#      Laboratoire d'Anthropologie et de Génétique des peuplements
#      Departement de génétique et d'évolution - UNiGE
# Date: 04 / 02 / 2016
# Language: Python 23
# 
# ---------------------------------------------------------------------------- #

# Import des modules nécéssaires
from fastaAndFurious import *
from Bio import pairwise2
import os

# La fonction qui calcul la distance de Hamming entre deux séquence
def dH(sequence_1, sequence_2, data_type='fasta_obj', force_align=1, align_type='system', gap=-1, missmatch=-1):
    '''
    Compute Hamming distance between 2 aligned sequence stored in a sequence class
    object.

    <data_type> is either (default) a fasta_obj class, or if 'nucl' nucleotidic
    sequences. 'fasta_obj' by default is to avoid problem with scripts using
    an oldest version of dH (which was by this time based solely on fasta obj)

    By using <force_align>, the user make the automatic alignement (by clustalo)
    automatic at each run. Disable (set to 0) if the sequences are already aligned
    to save time by reducing computational load.

        <align_type> is either 'system' (alignement done by an external tool) or
        'biopython' (in this case, alignement is done with biopython internal function).
        
        <gap> and <missmatch> are gap's and missmatch's costs.
    '''
    # Reconnaissance des types de données fournies
    if data_type  == 'nucl':
        name1 = 'SEQUENCE 1'
        seq1  = sequence_1
    
        name2 = 'SEQUENCE 2'
        seq2  = sequence_2

    if data_type  == 'fasta_obj':
        name1 = sequence_1.header[1:-1]
        seq1  = sequence_1.seq
    
        name2 = sequence_2.header[1:-1]
        seq2  = sequence_2.seq
    # Si l'utilisateur le demande, l'alignement est effectué
    if force_align == 1:
        if align_type == 'system':
            writeFasta([fasta_obj(name1, ungap(seq1)), fasta_obj(name2, ungap(seq2))], 'dH_TEMP.1')
            os.system('clustalo -i dH_TEMP.1 -o dH_TEMP.2 --force')
            data  = readFasta('dH_TEMP.2')
            
            name1 = data[0].header[1:-1]
            seq1  = data[0].seq
            
            name2 = data[1].header[1:-1]
            seq2  = data[1].seq
                        
            os.system('rm dH_TEMP.[1-2]')
                        
        elif align_type == 'biopython':
            bio_align = pairwise2.align.globalxs(ungap(seq1),ungap(seq2), missmatch, gap)[0]
            seq1      = bio_align[0] 
            seq2      = bio_align[1]

    if len(seq1) != len(seq2):
        print('Error: length differ\t%s\t<->\t%s' % (name1, name2))
        return

    # Computation of the Hamming Distance between the two sequences
    dist = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

# Fonction pour calculer les distances de Hamming entre plusieurs séquences
def multipleHamming(infile):
    '''
    Compute Hamming distances between all the sequences in a fasta file and write a
    distance matrix as output. Automatically align sequences with clustalo.
    
    Alignement is performed beforethe distance computation (saving time).
    '''
    data = readFasta(infile)
    
    data = alignSeqs(data)
    
    outfile = '%s_multipleHamming.txt' % infile.split('.')[0]
    outF = open(outfile, 'w')
    
    # headers
    for seq in data:
        outF.write('\t%s' % seq.header[1:])
    
    # calculs
    for seq1 in data:
        outF.write('\n%s' % seq1.header[1:])
        for seq2 in data:
            outF.write('\t%s' % dH(seq1, seq2, force_align = 0))

    outF.close()

# Retire les gaps dans une séquence
def ungap(sequence, filtered='-'):
    '''
    À partir d'une séquence fournie, renvoie ls même séquence sans les gap.
    
        filtered = Retirer aussi les N et |

    Version 2 - Plus de boucle, utilisation de methodes de la class <str>.
        Plus rapide.
    '''
    ret_seq = sequence
    
    for char in filtered:
        ret_seq = ''.join(ret_seq.split(char))
    
    return ret_seq

# Une fonction pour aligner l'ensemble des séquences d'une liste d'objets fasta_obj
def alignSeqs(liste_fastaObj, method='clustalo'):
    '''
    Aligne entre elles l'ensemble des séquences contenues dans une liste d'objets de
    type fasta_obj.
    method: Préciser l'algorithme d'alignement à utiliser. ClustalO est bien en multiple,
        Pour du pairwise ou des séquences similaires preferer muscle. 
    '''
    try:
        rm_output_1 = os.system('rm alignSeqs.1')
    except: 
        pass
    
    try:
        rm_output_2 = os.system('rm alignSeqs.2')
    except:
        pass
    
    # supression des gaps pour le bon fonctionnement de l'aligneur
    for item in liste_fastaObj:
        item.seq = ungap(item.seq)

    # operation d'I/O
    writeFasta(liste_fastaObj, 'alignSeqs.1')

    # alignement
    if method == 'clustalo':
        std_output = os.system('clustalo -i alignSeqs.1 -o alignSeqs.2 --force')
    elif method == 'muscle':
        std_output = os.system('muscle -in  alignSeqs.1 -out  alignSeqs.2 -quiet -stable')
    else:
        print('Methode Inconnue')
        return 0

    # verification du bon deroulement de l'alignement
    if std_output != 0:
        return std_output 

    # lecture des sequences alignées et supression des fichiers temporaires
    liste_aligned = readFasta('alignSeqs.2')
    std_output = os.system('rm alignSeqs.[1-2]')

    return liste_aligned

#EOF
