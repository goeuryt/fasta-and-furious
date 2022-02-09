#! /usr/bin/python3
# -*-coding:Utf-8 -*

# ---------------------------------------------------------------------------- #
#
# Les différentes fonction pour interagir avec un fichier fasta.
#
# ---------------------------------------------------------------------------- #
#
# Auteur: Thomas Goeury
#	  Laboratoire d'Anthropologie et de Génétique des peuplements
#	  Departement de génétique et d'évolution - UNiGE
# Date: 17 / 13 / 2016 
# Language: Python3
# 
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
#				   *					       
#				CLASSES					       
#				   *					       
# ---------------------------------------------------------------------------- #
class Fasta_obj:
    '''
    Class minimaliste pour stocker les infos d'une séquence fasta.
    
    Sert juste à fournir les infos pour un traitement ultérieur (identification
    d'infos supplémentaires dans les headers par ex.)
    
    Attributs:
      .header  -> Les infos dans les headers de la séquences
      .seq     -> La séquence nucléotidique
      .sup     -> Argument supplémentaire (si besoin est)
    Méthodes:
      .ungap() -> Retire les gaps de la séquence (ou autre caractères
                  si spécifiés)
        '''
    def __init__(self, info, nucleotids):
        self.header = info
        self.seq    = nucleotids
        self.sup    = 0
        
    def ungap(self, filtered = '-'):
        ''' Retire les gaps des séquences'''
        from fasTools import ungap as fasToolsUngap
        self.seq = fasToolsUngap(self.seq, filtered)


# ---------------------------------------------------------------------------- #
#				    *					       #
#				FONCTIONS				       #
#				    *					       #
# ---------------------------------------------------------------------------- #

#
## OUVERTURE DE FICHIERS FASTA
#
def openFile(fichier):
    '''
    A partir d'un chemin donné en argument, permet de récupérer l'information 
    contenue dedans sous formes d'une liste où chaque item est une ligne du fichier.
    '''
    FICHIER = open(fichier, 'r')
    DATA    = FICHIER.readlines()
    FICHIER.close()
    return DATA

def findSeq(fasta):
    '''
    Cette fonction recherche dans tout le texte la position des phrases commençant 
    par '>', qui correspondent donc aux lignes d'informations, puis à partir de là 
    calcul les positions de chaque séquence nucléotidique pour les extraires 
    '''
    # Recherche des positions de début de séquences
    positions = list()
    seq_list = list()
    for i in range(0, len(fasta)):
            if fasta[i].find('>') == 0:
                    positions.append(i)
    # Extraction des informations
    for i in range(0, len(positions)) :
        nucleotids = str()
        nucleotids_temp = str()
        infos = fasta[positions[i]][:-1] # Extraction des tags
        if i == len(positions) - 1: # Cas particulier où on en est au dernier élément de la liste
            nucleotids_temp = fasta[positions[i] + 1 : len(fasta)]
        else:
            nucleotids_temp = fasta[positions[i] + 1 : positions[i+1]]
        for j in nucleotids_temp:
            nucleotids = '%s%s' % (nucleotids, j[:-1])
        # Création d'un objet de class fasta pour stocker la séquence
        seq_temp = Fasta_obj(infos, nucleotids)
        # Ajout de cet objet à la liste des objets/sequences déjà trouvées
        seq_list.append(seq_temp)
    return seq_list

def readFasta(path):
    '''
    Permet d'utiliser en une seule fonction l'ensemble des fonctions au dessus
    '''
    return findSeq(openFile(path))

def fastaToDic(liste, key = 'header'):
    '''
    Transforme une liste d'objets "fasta_obj" en un dictionnaire.
    <Key> détermine si le header ou la sequence sera utilisée en tant
    que clef du dictionnaire. L'autre étant la valeur contenue dans 
    l'entrée.
    '''
    fastaDic = {}
    if key == 'header':
        for item in liste:
            try: 
                fastaDic[item.header].append(item.seq)
            except:
                fastaDic[item.header] = item.seq
    elif key == 'sequence':
        for  item in liste:
            try: 
                fastaDic[item.seq].append(item.header)
            except:
                fastaDic[item.seq] = item.header
    return fastaDic

#
## ECRITURE DE FICHIERS FASTA
#

def writeSequence(seq, fichier, endline):
    ''' 
    Ecrit correctement les séquences en suivant les normes fasta (séquences
    nucléotidiques sur 80 caractère max).
    '''
    from math import ceil
    
    nb_iterations=int(ceil(len(seq)/80.0))
    
    for i in range(0, nb_iterations):
        fichier.write('%s%s' % (seq[(i*80):(i+1)*80], endline))

# Ecrit un fichier sample sur le disque
def writeFasta(sample, file_name, mode='w'):
    ''' 
    Ecrit dans le fichier file_name l'ensemble des séquences contenue dans sample,
    sous forme FASTA.
    Sample doit être une liste d'objet de class <fasta_obj> tel que renvoyé par 
    readFasta().
    '''
    endline = '\n'
    out = open(file_name, mode)
    for i in range(0, len(sample)):
        if i == len(sample): 
                endline = ''
        out.write('>%s\n' % sample[i].header.split('>')[-1])
        writeSequence(sample[i].seq, out, endline)
    out.close()

# Permet de demander à l'utilisateur le chemin vers un fichier fasta (ou multi fasta)
def askForFasta():
    '''
    Met le programme en pause et demande un nom de fichier valide à l'utilisateur.
    Tant que le nom de fichier n'est pas correct ou que le fichier n'est pas lisible
    (par ex. ce n'est pas un fichier fasta), le programme va continuer de demander à
    l'utilisateur un nom de fichier valide.
    Renvoie ce fichier sous forme d'une liste d'objets de classe <sequence>.
    '''
    liste_alleles = []
    while len(liste_alleles) == 0:
        print("Fasta file path:")
        alleles_path = raw_input() # Permet de demander un input à l'utilisateur qui sera convertit automatiquement en <str>
        try:
            liste_alleles = readFasta(alleles_path)
        except:
            print("** Error: '%s' doesn't exist **" % alleles_path)
    return liste_alleles

# Nettoye un fichier fasta de tout ce qui n'est pas [ATCG] 	
def cleanFastaFile(PATH, filters = '-|N'):
    '''Nettoye le fichier fasta de tout ce qui est passé en argument de "filtered" '''
    fasta  = readFasta(PATH)
    for f in fasta:
        f.ungap(filters)
    writeFasta(fasta, PATH)





# EOF
