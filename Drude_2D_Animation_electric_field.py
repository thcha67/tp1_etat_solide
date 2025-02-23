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
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--efield', type=float, default=0.1, help='Electric field in V/M')
parser.add_argument("--tau", type=float, default=200e-6, help="Time constant for the slowing down based on the Maxwell distribution of speeds")
parser.add_argument("--noshow", action="store_true", default=False, help="Don't show the animation")
args = parser.parse_args()


# win = 500 # peut aider à définir la taille d'un autre objet visuel comme un histogramme proportionnellement à la taille du canevas.

# Déclaration de variables influençant le temps d'exécution de la simulation
Nelectrons = 500  # change this to have more or fewer electrons
dt = 1E-7  # pas d'incrémentation temporel
Nions = 16  # nombre de coeurs ioniques

# Déclaration de variables physiques "Typical values"
DIM = 2 #Nombre de degrés de liberté de la simulation 
#mass = 4E-3/6E23 # helium mass
mass = 9.10938356*10**(-31) # masse de l'électron
Relectron = 0.01 # wildly exaggerated size of an electron
Rion = 0.01 # wildly exaggerated size of an atom
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature

_tau = args.tau # constante de temps du au ralentissement basé sur la distribution de Maxwell des vitesses
q = 1.60217662*10**(-19)  # charge élémentaire
electric_field = args.efield  # champ électrique

#### CANEVAS DE FOND ####
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container and spheres below
if not args.noshow:
    animation = canvas(width=750, height=500) # , align='left')
    animation.range = L
# animation.title = 'Théorie cinétique des gaz parfaits'
# s = """  Simulation de particules modélisées en sphères dures pour représenter leur trajectoire ballistique avec collisions. Une sphère est colorée et grossie seulement pour l’effet visuel permettant de suivre sa trajectoire plus facilement dans l'animation, sa cinétique est identique à toutes les autres particules.

# """
# animation.caption = s

#### ARÊTES DE BOÎTE 2D ####
d = L/2+Relectron
r = 0.005
cadre = curve(color=gray, radius=r)
cadre.append([vector(-d,-d,0), vector(d,-d,0), vector(d,d,0), vector(-d,d,0), vector(-d,-d,0)])

#### POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES ####
Pink_Electron_Index = 0 # garde une sphère plus grosse et colorée parmis toutes les grises
List_Electrons_Obj: list[simple_sphere] = [] # Objet qui contiendra les sphères pour l'animation
List_Momentum_Electrons: list[vector] = [] # quantité de mouvement des sphères
List_Position_Electrons: list[vector] = [] # position des sphères
List_Index_Electrons: list[int] = [] # index des sphères
pavg = sqrt(2*mass*(DIM/2)*k*T) # average kinetic energy in 3D p**2/(2mass) = (3/2)kT : Principe de l'équipartition de l'énergie en thermodynamique statistique classique


for electron in range(Nelectrons):
    # liste de l'index de chaque sphère
    List_Index_Electrons.append(electron)

    # position aléatoire qui tient compte que l'origine est au centre de la boîte
    x = L*random()-L/2 
    y = L*random()-L/2
    z = 0
    if electron == Pink_Electron_Index:  # garde une sphère plus grosse et colorée parmis toutes les grises
        List_Electrons_Obj.append(simple_sphere(pos=vector(x,y,z), radius=0.03, color=color.magenta)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: 
        List_Electrons_Obj.append(simple_sphere(pos=vector(x,y,z), radius=Relectron, color=gray))


    # liste de la position initiale de toutes les sphères
    List_Position_Electrons.append(vec(x,y,z))

    # theta = pi*random() # direction de coordonnées sphériques, superflue en 2D
    phi = 2*pi*random() # direction aléatoire pour la quantité de mouvement
    px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
    py = pavg*sin(phi)
    pz = 0
    List_Momentum_Electrons.append(vector(px,py,pz)) # liste de la quantité de mouvement initiale de toutes les sphères


for ion in range(Nions):
    idx = electron + ion
    # reparti les coeurs ioniques dans la boîte uniformément
    assert np.sqrt(Nions) % 1 == 0, "Le nombre de coeurs ioniques doit être un carré parfait"
    x_pos_ion = (L/2 - (ion % np.sqrt(Nions)) * L/np.sqrt(Nions)) - ((L/2) / np.sqrt(Nions))
    y_pos_ion = (L/2 - (ion // np.sqrt(Nions)) * L/np.sqrt(Nions)) - ((L/2) / np.sqrt(Nions))
    z_pos_ion = 0
    List_Index_Electrons.append(idx)
    List_Electrons_Obj.append(simple_sphere(pos=vector(x_pos_ion,y_pos_ion,z_pos_ion), radius=0.03, color=color.red)) #, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    List_Position_Electrons.append(vec(x_pos_ion,y_pos_ion,z_pos_ion))
    List_Momentum_Electrons.append(vector(0,0,0))

#### FONCTION POUR IDENTIFIER LES COLLISIONS, I.E. LORSQUE LA DISTANCE ENTRE LES CENTRES DE 2 SPHÈRES EST À LA LIMITE DE S'INTERPÉNÉTRER ####
def checkCollisions() -> list[tuple[int, int]]:
    electrons_with_hits = []
    ions_with_hits = []
    r2 = Relectron + Rion   # distance critique où les 2 sphères entre en contact à la limite de leur rayon
    r2 *= r2   # produit scalaire pour éviter une comparaison vectorielle ci-dessous
    for i in range(Nelectrons):
        ai = List_Position_Electrons[i]
        for j in range(Nions):  #On compare chaque atome avec chaque coeur ionique
            aj = List_Position_Electrons[j + Nelectrons]
            dr = ai - aj
            if mag2(dr) < r2:
                electrons_with_hits.append(i)
                ions_with_hits.append(j + Nelectrons)

    hitlist = zip(electrons_with_hits, ions_with_hits)
    no_hitlist = list(set(range(Nelectrons)) - set(electrons_with_hits))

    return hitlist, no_hitlist

#### BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt ####
## ATTENTION : la boucle laisse aller l'animation aussi longtemps que souhaité, assurez-vous de savoir comment interrompre vous-même correctement (souvent `ctrl+c`, mais peut varier)
## ALTERNATIVE : vous pouvez bien sûr remplacer la boucle "while" par une boucle "for" avec un nombre d'itérations suffisant pour obtenir une bonne distribution statistique à l'équilibre

momentum_averages = [] # quantité de mouvement moyenne (norme)
position_x_averages = [] # position moyenne en x
position_y_averages = [] # position moyenne en y

j = 10000 # nombre d'itérations
for _ in tqdm(range(j)):
    #rate(300)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

    #### DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
    List_Vitesse_Atoms: list[vector] = []   # vitesse instantanée de chaque sphère
    List_Ajout_Distance: list[vector] = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
    for electron in List_Index_Electrons:
        if electron < Nelectrons:
            List_Vitesse_Atoms.append(List_Momentum_Electrons[electron]/mass)   # par définition de la quantité de nouvement pour chaque sphère
            List_Ajout_Distance.append(List_Vitesse_Atoms[electron] * dt)   # différence avant pour calculer l'incrément de position

            # nouvelle position de l'atome après l'incrément de temps dt
            List_Electrons_Obj[electron].pos = List_Position_Electrons[electron] = List_Position_Electrons[electron] + List_Ajout_Distance[electron] 
        else:
            List_Vitesse_Atoms.append(vector(0,0,0))

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE ####
    for electron in List_Index_Electrons:
        loc = List_Position_Electrons[electron]
        if loc.x > L/2 and List_Momentum_Electrons[electron].x > 0:
            List_Momentum_Electrons[electron].x *= -1
        if loc.x < -L/2 and List_Momentum_Electrons[electron].x < 0:
            List_Momentum_Electrons[electron].x *= -1

        if loc.y > L/2 and List_Momentum_Electrons[electron].y > 0:
            List_Momentum_Electrons[electron].y *= -1
        if loc.y < -L/2 and List_Momentum_Electrons[electron].y < 0:
            List_Momentum_Electrons[electron].y *= -1

    momentum_average = 0
    for idx in range(Nelectrons):
        momentum_average += mag(List_Momentum_Electrons[idx])
    momentum_average /= Nelectrons

    T = (momentum_average**2) / (mass * DIM * k)  # température moyenne de l'ensemble de particules

    #### LET'S FIND THESE COLLISIONS!!! ####
    hitlist, no_hitlist = checkCollisions() # list[tuple[int, int]], set[int]

    #### CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES ####
    for hit in hitlist:
        # définition de nouvelles variables pour chaque paire de sphères en collision
        # extraction du numéro des 2 sphères impliquées à cette itération
        electron = hit[0]
        coeur = hit[1]

        pos_electron = List_Position_Electrons[electron]
        pos_coeur = List_Position_Electrons[coeur]
        momentum_electron = List_Momentum_Electrons[electron]
        momentum_electron_norm = mag(momentum_electron)

        new_momentum_electron_norm = mass * np.random.rayleigh(scale=np.sqrt(k * T / mass), size=1)  # collision inélastique avec le coeur ionique

        theta = np.random.uniform(0, 2 * np.pi)
        new_p_x_electron = new_momentum_electron_norm * np.cos(theta)
        new_p_y_electron = new_momentum_electron_norm * np.sin(theta)
        new_p_z_electron = 0

        new_momentum_electron = vector(new_p_x_electron.item(), new_p_y_electron.item(), new_p_z_electron)

        List_Momentum_Electrons[electron] = new_momentum_electron

        List_Momentum_Electrons[coeur] = vector(0,0,0)  # le coeur ionique reste immobile
        
        new_position = pos_electron + List_Momentum_Electrons[electron]/mass * dt

        List_Position_Electrons[electron] = new_position

    for no_hit in no_hitlist:  
        momentum_electron = List_Momentum_Electrons[no_hit] # get the momentum of the electron
        momentum_electron_x = momentum_electron.x # get the x component of the momentum
        momentum_electron_y = momentum_electron.y
        momentum_electron_z = momentum_electron.z

        momentum_increment = dt*(-momentum_electron_x/_tau + q*electric_field) # apply Coulomb's force as in A/M 1.12

        position = List_Position_Electrons[no_hit]

        if position.x > L/2:
            if momentum_increment > 0:
                momentum_electron_x = 0
                List_Position_Electrons[no_hit].x = L/2
        elif position.x < -L/2:
            if momentum_increment < 0:
                momentum_electron_x = 0
                List_Position_Electrons[no_hit].x = -L/2
        else:
            momentum_electron_x += momentum_increment

            List_Momentum_Electrons[no_hit] = vector(momentum_electron_x, momentum_electron_y, momentum_electron_z)

            new_position = position + List_Momentum_Electrons[no_hit]/mass * dt

            List_Position_Electrons[no_hit] = new_position


    momentum_averages.append(momentum_average)
    position_x_averages.append(np.mean([List_Position_Electrons[i].x for i in range(Nelectrons)]))
    position_y_averages.append(np.mean([List_Position_Electrons[i].y for i in range(Nelectrons)]))
