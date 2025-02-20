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
from scipy.stats import maxwell

# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Natoms = 200  # change this to have more or fewer atoms
dt = 1E-7  # pas d'incrémentation temporel
Ncoeurs = 16  # nombre de coeurs ioniques

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
#mass = 4E-3/6E23 # helium mass
mass = 9.10938356*10**(-31) # masse de l'électron
Ratom = 0.01 # wildly exaggerated size of an atom
Rcoeur = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
animation = canvas(width=750, height=500) # , align='left')
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


for coeur_1 in range(Ncoeurs):
    idx = Atom_1 + coeur_1
    # reparti les coeurs ioniques dans la boîte uniformément
    assert np.sqrt(Ncoeurs) % 1 == 0, "Le nombre de coeurs ioniques doit être un carré parfait"
    x_ionique = (L/2 - (coeur_1 % np.sqrt(Ncoeurs)) * L/np.sqrt(Ncoeurs)) - ((L/2) / np.sqrt(Ncoeurs))
    y_ionique = (L/2 - (coeur_1 // np.sqrt(Ncoeurs)) * L/np.sqrt(Ncoeurs)) - ((L/2) / np.sqrt(Ncoeurs))
    z_ionique = 0
    List_Index_Atom.append(idx)
    List_Atoms_Obj.append(simple_sphere(pos=vector(x_ionique,y_ionique,z_ionique), radius=0.03, color=color.red)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    List_Position_Atom.append(vec(x_ionique,y_ionique,z_ionique))
    List_Momentum_Atoms.append(vector(0,0,0))

#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions() -> list[tuple[int, int]]:
    hitlist = []   # initialisation
    r2 = Ratom + Rcoeur   # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Natoms):
        ai = List_Position_Atom[i]
        for j in range(Ncoeurs):  #On compare chaque atome avec chaque coeur ionique
            aj = List_Position_Atom[j + Natoms]
            dr = ai - aj
            if mag2(dr) < r2:
                hitlist.append((i,j + Natoms))

    return hitlist

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

c = 0 
while c < 2000:
    rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    List_Vitesse_Atoms: list[vector] = []   # vitesse instantanée de chaque sphère
    List_Ajout_Distance: list[vector] = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for Atom_1 in List_Index_Atom:
        if Atom_1 < Natoms:
            List_Vitesse_Atoms.append(List_Momentum_Atoms[Atom_1]/mass)   # par définition de la quantité de nouvement pour chaque sphère
            List_Ajout_Distance.append(List_Vitesse_Atoms[Atom_1] * dt)   # différence avant pour calculer l'incrément de position

            # nouvelle position de l'atome après l'incrément de temps dt
            List_Atoms_Obj[Atom_1].pos = List_Position_Atom[Atom_1] = List_Position_Atom[Atom_1] + List_Ajout_Distance[Atom_1] 
        else:
            List_Vitesse_Atoms.append(vector(0,0,0))



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


    momentum_average_squared = 0
    for idx in range(Natoms):
        momentum_average_squared += mag(List_Momentum_Atoms[idx])
    momentum_average_squared /= Natoms
    T = (momentum_average_squared**2) / (mass * DIM * k)  # température moyenne de l'ensemble de particules
    print(T)

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist: list[tuple[int, int]] = checkCollisions()


    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for hit in hitlist:

        # définition de nouvelles variables pour chaque paire de sphères en collision
        # extraction du numéro des 2 sphères impliquées à cette itération
        electron = hit[0]  
        coeur = hit[1]

        pos_electron = List_Position_Atom[electron]
        pos_coeur = List_Position_Atom[coeur]
        momentum_electron = List_Momentum_Atoms[electron]
        momentum_electron_norm = mag(momentum_electron)

        new_momentum_electron_norm = mass * np.random.rayleigh(scale=np.sqrt(k * T / mass), size=1)  # collision inélastique avec le coeur ionique

        theta = np.random.uniform(0, 2 * np.pi)
        new_p_x_electron = new_momentum_electron_norm * np.cos(theta)
        new_p_y_electron = new_momentum_electron_norm * np.sin(theta)
        new_p_z_electron = 0
        new_momentum_electron = vector(new_p_x_electron, new_p_y_electron, new_p_z_electron)

        List_Momentum_Atoms[electron] = new_momentum_electron
        List_Momentum_Atoms[coeur] = vector(0,0,0)  # le coeur ionique reste immobile

        List_Position_Atom[electron] = pos_electron + List_Momentum_Atoms[electron]/mass * dt


    c += 1
