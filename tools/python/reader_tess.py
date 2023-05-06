# -*- coding:utf8 -*-
#
# Reader for neper outputs
# Author: L. Moutin
# Date: 25/10/2021
#
# Copyright : see License.txt


# ------------------- Fonctions secondaires ---------------------------------
# Pour la lecture des coordonées du domaine
def lecture_coord(dim,ligne,sortie):
    # Cas 2D
    if (dim == 2):
        if ligne.find('x1y0')>0: sortie.write("'Lx'"+" "+ligne[6:20]+"\n")
        if ligne.find('x0y1')>0: sortie.write("'Ly'"+" "+ligne[21:35]+"\n")
    # Cas 3D
    if (dim == 3):
        if ligne.find('x1y0z0')>0: sortie.write("'Lx'"+" "+ligne[6:20]+"\n")
        if ligne.find('x0y1z0')>0: sortie.write("'Ly'"+" "+ligne[21:35]+"\n")
        if ligne.find('x0y0z1')>0: sortie.write("'Lz'"+" "+ligne[36:50]+"\n")
    return




# ------------------- Fonction principale ---------------------------------
def reader_tess(FichierEntree):

    FichierEntree = './'+FichierEntree+'.tess'
    FichierSortie = 'FichierSeeds'+'.txt'

    # Ouverture des fichiers d'entrée et de sortie
    with open(FichierEntree,'r') as entree:
        with open(FichierSortie,'w') as sortie:
            lignes_entree = entree.readlines()
            i = 0; dim = 0
            # Premier parcours pour lire les infos : Dimension, nombre de cellues et coordonnées du domaine
            for ligne in lignes_entree:
                # lecture de la dimension
                if ligne.find('**general')>0:
                    sortie.write("'Dimension'"+" "+lignes_entree[i+1][3:4]+"\n")
                    dim = int(lignes_entree[i+1][3:4])
               # lecture du nombre de cellules
                if ligne.find('*cell')>0:
                    sortie.write("'Points'"+" "+lignes_entree[i+1])
            # lecture de la taille du domaine (Lx, Ly, Lz)
                lecture_coord(dim,ligne,sortie)

                i += 1



# Second parcours pour récuperer les infos sur les seeds
            cpt_s = False
            for ligne in lignes_entree:
                if (cpt_s == True):
                    if ligne.find('*')>0:
                        return
                    sortie.write(ligne[5:])
                if ligne.find('*seed')>0:
                    cpt_s = True

    return

