#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 200  # change this to have more or fewer atoms
dt = 1E-5  # pas d'incrémentation temporel

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
mass = 4E-3/6E23 # helium mass
Ratom = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas( width=750, height=500) # , align='left')
animation.range = L
# animation.title = 'Théorie cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Ratom
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Pink_Atom_Index = 0 # garde une sphère plus grosse et colorée parmis toutes les grises
List_Atoms_Obj: list[simple_sphere] = [] # Objet qui contiendra les sphères pour l'animation
List_Momentum_Atoms: list[vector] = [] # quantité de mouvement des sphères
List_Position_Atom: list[vector] = [] # position des sphères
List_Index_Atom: list[int] = [] # index des sphères
pavg = sqrt(2*mass*(DIM/2)*k*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique

for Atom_1 in range(Natoms):
    # liste de l'index de chaque sphère
    List_Index_Atom.append(Atom_1)

    # position aléatoire qui tient compte que l'origine est au centre de la boîte
    x = L*random()-L/2 
    y = L*random()-L/2
    z = 0
    if Atom_1 == Pink_Atom_Index:  # garde une sphère plus grosse et colorée parmis toutes les grises
        List_Atoms_Obj.append(simple_sphere(pos=vector(x,y,z), radius=0.03, color=color.magenta)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: 
        List_Atoms_Obj.append(simple_sphere(pos=vector(x,y,z), radius=Ratom, color=gray))

    # liste de la position initiale de toutes les sphères
    List_Position_Atom.append(vec(x,y,z))

    # theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    List_Momentum_Atoms.append(vector(px,py,pz)) # liste de la quantité de mouvement initiale de toutes les sphères

#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions() -> list[tuple[int, int]]:
    hitlist = []   # initialisation
    r2 = 2*Ratom   # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        ai = List_Position_Atom[i]
        for j in range(i) :
            aj = List_Position_Atom[j]
            dr = ai - aj   # la boucle dans une boucle itère pour calculer cette distance vectorielle dr entre chaque paire de sphère
            if mag2(dr) < r2:   # test de collision où mag2(dr) qui retourne la norme élevée au carré de la distance intersphère dr
                hitlist.append((i,j)) # liste numérotant toutes les paires de sphères en collision
    return hitlist

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

c = 0 
while c < 200:
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    List_Vitesse_Atoms: list[vector] = []   # vitesse instantanée de chaque sphère
    List_Ajout_Distance: list[vector] = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for Atom_1 in List_Index_Atom:
        List_Vitesse_Atoms.append(List_Momentum_Atoms[Atom_1]/mass)   # par définition de la quantité de nouvement pour chaque sphère
        List_Ajout_Distance.append(List_Vitesse_Atoms[Atom_1] * dt)   # différence avant pour calculer l'incrément de position

        # nouvelle position de l'atome après l'incrément de temps dt
        List_Atoms_Obj[Atom_1].pos = List_Position_Atom[Atom_1] = List_Position_Atom[Atom_1] + List_Ajout_Distance[Atom_1]  

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for Atom_1 in List_Index_Atom:
        loc = List_Position_Atom[Atom_1]

        if abs(loc.x) > L/2:
            if loc.x < 0:
                List_Momentum_Atoms[Atom_1].x = abs(List_Momentum_Atoms[Atom_1].x)  # renverse composante x au mur de gauche
            else: 
                List_Momentum_Atoms[Atom_1].x = -abs(List_Momentum_Atoms[Atom_1].x)   # renverse composante x au mur de droite

        if abs(loc.y) > L/2:
            if loc.y < 0: 
                List_Momentum_Atoms[Atom_1].y = abs(List_Momentum_Atoms[Atom_1].y)  # renverse composante y au mur du bas
            else: 
                List_Momentum_Atoms[Atom_1].y =  -abs(List_Momentum_Atoms[Atom_1].y)  # renverse composante y au mur du haut


    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist: list[tuple[int, int]] = checkCollisions()


    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for hit in hitlist:

        # définition de nouvelles variables pour chaque paire de sphères en collision
        # extraction du numéro des 2 sphères impliquées à cette itération
        Atom_1 = hit[0]  
        Atom_2 = hit[1]

        ptot = List_Momentum_Atoms[Atom_1]+List_Momentum_Atoms[Atom_2]   # quantité de mouvement totale des 2 sphères
        mtot = 2*mass    # masse totale des 2 sphères
        Vcom = ptot/mtot   # vitesse du référentiel barycentrique/center-of-momentum (com) frame
        posi = List_Position_Atom[Atom_1]   # position de chacune des 2 sphères
        posj = List_Position_Atom[Atom_2]
        vi = List_Momentum_Atoms[Atom_1]/mass   # vitesse de chacune des 2 sphères
        vj = List_Momentum_Atoms[Atom_2]/mass
        rrel = posi-posj  # vecteur pour la distance entre les centres des 2 sphères
        vrel = vj-vi   # vecteur pour la différence de vitesse entre les 2 sphères

        # exclusion de cas où il n'y a pas de changements à faire
        if vrel.mag2 == 0: continue  # exactly same velocities si et seulement si le vecteur vrel devient nul, la trajectoire des 2 sphères continue alors côte à côte
        if rrel.mag > Ratom: continue  # one atom went all the way through another, la collision a été "manquée" à l'intérieur du pas deltax

        # calcule la distance et temps d'interpénétration des sphères dures qui ne doit pas se produire dans ce modèle
        dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta) où theta is the angle between vrel and rrel:
        dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        alpha = asin(dy/(2*Ratom))  # alpha is the angle of the triangle composed of rrel, path of atom j, and a line from the center of atom i to the center of atom j where atome j hits atom i
        d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        deltat = d/vrel.mag         # time spent moving from first contact to position inside atom

        #### CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION ####
        posi = posi-vi*deltat   # back up to contact configuration
        posj = posj-vj*deltat
        pcomi = List_Momentum_Atoms[Atom_1]-mass*Vcom  # transform momenta to center-of-momentum (com) frame
        pcomj = List_Momentum_Atoms[Atom_2]-mass*Vcom
        rrel = hat(rrel)    # vecteur unitaire aligné avec rrel
        pcomi = pcomi-2*dot(pcomi,rrel)*rrel # bounce in center-of-momentum (com) frame
        pcomj = pcomj-2*dot(pcomj,rrel)*rrel
        List_Momentum_Atoms[Atom_1] = pcomi+mass*Vcom # transform momenta back to lab frame
        List_Momentum_Atoms[Atom_2] = pcomj+mass*Vcom
        List_Position_Atom[Atom_1] = posi+(List_Momentum_Atoms[Atom_1]/mass)*deltat # move forward deltat in time, ramenant au même temps où sont rendues les autres sphères dans l'itération
        List_Position_Atom[Atom_2] = posj+(List_Momentum_Atoms[Atom_2]/mass)*deltat

    c += 1
