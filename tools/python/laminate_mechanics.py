# -*- coding:utf8 -*-
#
# Modeles homogenes de lamines mecaniques
# Voir Correcteurs_lamines.pdf
# Date: 10/11/2020
#
# Copyright : see License.txt
#
# Descriptif:
# On a un lamine de matériaux isotropes selon l'axe des x
# | mat 1 | mat 2 | mat 3 | ... | mat n |
# | mat 1 | mat 2 | mat 3 | ... | mat n |
# | mat 1 | mat 2 | mat 3 | ... | mat n |
# | mat 1 | mat 2 | mat 3 | ... | mat n |
# ----------------------------------------> x
# On utilise la convention de Voigt
# (s1, s2, s3, s3, s4, s5, s6)  = C (eps1, eps2, eps3, 2 eps4, 2 eps5, 2 eps6).
# Voir https://www-git-cad.intra.cea.fr/DEC/collaboratif/jv253511/libpy
# 

import numpy as np
from math import floor


# format de l'entree
listeExemple = [{'lambda': 100., 'mu': 75., 'fracVol': 0.5}, {'lambda': 200., 'mu': 150., 'fracVol': 0.5}]
# lambda, mu  -> constantes mecaniques
# fracVol -> fraction volumique d'une phase


# BEGIN{Fonctions auxiliaires}
def moyenne(liste, fonction):
    return sum([mat['fracVol'] * fonction(mat) for mat in liste])


def fonction1(mat):
    return 1. / (mat['lambda'] + 2 * mat['mu'])


def fonction2(mat):
    return mat['lambda'] / (mat['lambda'] + 2 * mat['mu'])


def fonction3(mat):
    return 1. / mat['mu']
# END{Fonctions auxiliaires}


# BEGIN{Parties non-triviales des deformations issues de deplacements imposes}
def epsilon1xx(liste, phase):
    return fonction1(liste[phase]) / moyenne(liste, fonction1)


def epsilon2xx(liste, phase):
    return moyenne(liste, fonction2) * fonction1(liste[phase]) / moyenne(liste, fonction1) - fonction2(liste[phase])


def epsilon3xx(liste, phase):
    return epsilon2xx(liste, phase)


def epsilon5xz(liste, phase):
    return fonction3(liste[phase]) / moyenne(liste, fonction3)


def epsilon6xy(liste, phase):
    return epsilon5xz(liste, phase)
# END{Parties non-triviales des deformations issues de deplacements imposes}


# Matrice de passage de déformation moyenne à déformation local en notations de Voigt
def matriceDePassage(liste, phase):
    return [[epsilon1xx(liste, phase), epsilon2xx(liste, phase), epsilon3xx(liste, phase), 0, 0, 0],
            [0, 1., 0, 0, 0, 0],
            [0, 0, 1., 0, 0, 0],
            [0, 0, 0, 1., 0, 0],
            [0, 0, 0, 0, epsilon5xz(liste, phase), 0],
            [0, 0, 0, 0, 0, epsilon6xy(liste, phase)]
            ]


# Matrice de rigidite en notations de Voigt
def matRigidite(liste, phase):
    mat = liste[phase]
    c1 = 2 * mat['mu'] + mat['lambda']
    c2 = mat['lambda']
    c3 = mat['mu']
    return [[c1, c2, c2, 0, 0, 0],
            [c2, c1, c2, 0, 0, 0],
            [c2, c2, c1, 0, 0, 0],
            [0, 0, 0, c3, 0, 0],
            [0, 0, 0, 0, c3, 0],
            [0, 0, 0, 0, 0, c3]]


# Matrice de rigidite par phase
def matRigiditeEq(liste, phase):
    matRig = np.array(matRigidite(liste, phase))
    matPass = np.array(matriceDePassage(liste, phase))
    res = np.dot(matRig, matPass)
    return list(res)


# Matrice homogeneisee
def matHomog(liste):
    matSomme = np.array([[0 for i in range(0, 6)] for j in range(0, 6)])
    for i in range(0, len(liste)):
        matSomme = matSomme + liste[i]['fracVol'] * np.array(matRigiditeEq(liste, i))
    return list(matSomme)


# Test de validation
def printMatHomog(liste):
    """
    Cet exemple est soutenu par le document .pdf cite en intro et un calcul via VER-TMFFT
    @Example
    >>> printMatHomog(listeExemple)
    [333, 133, 133, 0, 0, 0]
    [133, 368, 143, 0, 0, 0]
    [133, 143, 368, 0, 0, 0]
    [0, 0, 0, 112, 0, 0]
    [0, 0, 0, 0, 100, 0]
    [0, 0, 0, 0, 0, 100]
    """
    matHom = matHomog(liste)
    matHomAff = [[int(floor(elem)) for elem in ligne] for ligne in matHom]
    for ligne in matHomAff:
        print(ligne)


# printMatHomog(listeExemple)

if __name__ == '__main__':
    from doctest import testmod

    testmod(verbose=True)
