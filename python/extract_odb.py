# -*- coding: utf-8 -*-
"""
Extraction de la réponse harmonique depuis l'ODB du job harmonique
  - S11 (module complexe) au nœud RACINE, moyenné sur les éléments voisins
  - U3  (module complexe) au nœud BOUT
Sortie : python/reponse_1g.csv  (colonnes : freq_Hz, S11_MPa, U3_mm)

Lancement (au choix) :
  - Abaqus/CAE : File > Run Script... > choisir ce fichier
    (le répertoire de travail doit être vibro_fatigue/abaqus)
  - Abaqus Command : cd vers vibro_fatigue/abaqus puis
        abaqus python ../python/extract_odb.py
"""
from __future__ import print_function
import math
import os
from odbAccess import openOdb
from abaqusConstants import ELEMENT_NODAL

ODB_PATH = 'harmonique_2.odb'   # nom du job + .odb
CSV_PATH = os.path.join('..', 'python', 'reponse_1g.csv')
SET_RACINE = 'RACINE'
SET_BOUT = 'BOUT'


def trouver_set(odb, nom):
    """Cherche un node set au niveau assemblage, puis au niveau instance."""
    ra = odb.rootAssembly
    if nom in ra.nodeSets.keys():
        return ra.nodeSets[nom]
    for inst_nom in ra.instances.keys():
        inst = ra.instances[inst_nom]
        if nom in inst.nodeSets.keys():
            return inst.nodeSets[nom]
    raise KeyError("Node set '%s' introuvable dans l'ODB" % nom)


def module(re, im):
    return math.sqrt(re * re + im * im)


def composante_complexe(valeurs, idx):
    """Moyenne (partie réelle, partie imaginaire) d'une composante sur plusieurs valeurs."""
    # Boucle explicite : dans Abaqus/CAE, sum() est remplacé par la fonction
    # Abaqus de somme de champs, qui n'accepte pas les générateurs Python.
    re, im = 0.0, 0.0
    for v in valeurs:
        re += v.data[idx]
        if v.conjugateData is not None:
            im += v.conjugateData[idx]
    n = float(len(valeurs))
    return re / n, im / n


def main():
    odb = openOdb(ODB_PATH, readOnly=True)
    nom_step = list(odb.steps.keys())[-1]          # dernier step = harmonique
    step = odb.steps[nom_step]
    racine = trouver_set(odb, SET_RACINE)
    bout = trouver_set(odb, SET_BOUT)
    print("Step lu : %s  (%d frames)" % (nom_step, len(step.frames)))

    lignes = []
    for frame in step.frames:
        f = frame.frameValue                        # fréquence en Hz pour un step SSD
        if f <= 0.0:
            continue
        s_vals = frame.fieldOutputs['S'].getSubset(
            region=racine, position=ELEMENT_NODAL).values
        u_vals = frame.fieldOutputs['U'].getSubset(region=bout).values
        s11 = module(*composante_complexe(s_vals, 0))   # S11 = composante 0
        u3 = module(*composante_complexe(u_vals, 2))    # U3  = composante 2
        lignes.append((f, s11, u3))

    odb.close()
    lignes.sort()
    with open(CSV_PATH, 'w') as fic:
        fic.write('freq_Hz,S11_MPa,U3_mm\n')
        for f, s, u in lignes:
            fic.write('%.4f,%.6f,%.6f\n' % (f, s, u))

    f_pic, s_pic, _ = max(lignes, key=lambda l: l[1])
    print("%d points écrits dans %s" % (len(lignes), os.path.abspath(CSV_PATH)))
    print("Pic S11 : %.1f MPa à %.2f Hz" % (s_pic, f_pic))
    print("Pic U3  : %.2f mm" % max(l[2] for l in lignes))


main()
