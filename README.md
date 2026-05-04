# Routage

## Description

Ce projet consiste en la création d'un algorithme de routage météorologique pour les voiliers. Étant donné des condition météorologiques (champs de vent, courant) et les performances d'un bateau (fonction polaire), l'objectif est de trouver la route de plus court temps entre un point départ et un point d'arrivée. L'algorithme developpé est basé sur la méthode des isochrones.

![Screenshot](images/screenshot.png)


## Détail fichiers

vents.py : modèles de champs de vent V : x, y, t -> wind_dir, wind_speed

vents_grib.py : interpolations de prévisions de vent V : x, y, t -> wind_dir, wind_speed

polaires.py : polaires de vitesse de bateaux P : ang, wind_speed -> boat_speed

enveloppe.py : calcul de l'enveloppe d'un nuage de points

isochrone.py : mise en oeuvre de la méthode des isochrones, utilise un enveloppe.py

exec_.ipynb : exection et affichage de l'ensemble


## Roadmap

[x] Réduction temps de calcul par passage numpy, KDTree <br>
[x] Réduciotn temps de calcul par calcul partiel d'isochrones  <br>
[x] Interface avec barre de temps <br>
[x] Gestion côte fixe <br>
[x] Revoir pipeline donnée grib <br>
[x] Revoir pipeline donnée polaires <br>
[ ] Essayer nouv methode enveloppe k-plus proches voisins <br>

[ ] Passer les angles en rad <br>
[ ] Intégration courants <br>
[ ] Interpolation grib <br>
[ ] Fonction rafineur de route <br>
[ ] Corriger conversion def <-> nm <br>
[ ] Intégration marrée, côté non fixe <br>


## Questions & idées

- Quelle est la bonne transformation pour passer de distance à anlges latitude / longitude ? <br>
- Comment mieux gérer le points d'arrivée ? <br>