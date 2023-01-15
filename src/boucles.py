import os
import re
import argparse
from Bio.SeqUtils import seq1


parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store', dest='fichier_cible',
                    help='Nom du fichier ".out"', type=str, required=True)
parser.add_argument('-c', action='store', dest='chemin',
                    help='Chemin vers les sites modifiés dans le fichier ".out"', type=str, default="../Patch_Search/Sites_Modifies/")
parser.add_argument('-s', action='store', dest='score_seuil',
                    help='Seuil de score', type=float, default="16")
parser.add_argument('-l', action='store', dest='longueur_min',
                    help="Longueur minimale d'une boucle", type=int, default="4")    
parser.add_argument('-g', action='store', dest='taille_gap',
                    help='Seuil de tolérance aux gaps dans une boucle', type=int, default="0")


target = parser.parse_args().fichier_cible
chemin = parser.parse_args().chemin
score_seuil = float(parser.parse_args().score_seuil)
longueur_min = int(parser.parse_args().longueur_min)
gap = int(parser.parse_args().taille_gap)


if chemin[-1] != "/":
    chemin += "/"

if target[-4:] == ".out":
    target = target[:-4]

# Suppression des fichiers de sortie si ils existent déjà
if os.path.isfile(f"../res/{target}_boucles.fasta"):
    os.remove(f"../res/{target}_boucles.fasta")

if os.path.isfile(f"../res/scores.csv"):
    os.remove(f"../res/scores.csv")


# Pour chaque site de liasion avec un score supérieur au seuil défini, on récupère la liste de ses résidus
def Recuperation_residu(cible):

    # On va aussi créer un fichier csv avec tous les scores renvoyés
    with open(f"../res/{cible}.out", 'r') as target, open("../res/scores.csv", 'w') as data_score:
        data_score.write("Score\n")
        tab_score=[]
        for line in target:
            tableau_Residu = []

            # Permet d'enlever la première ligne ainsi que les rares lignes mal formatées
            if len(line) == 1 or len(line) < len(chemin) or 180 < len(line):
                line = "bad"
            if line.strip()[0] != "b":
                score=line.split()[-1]

                # On récupère la liste des scores où une similarité a été trouvée
                if float(score)!=0:
                    data_score.write(f"{score}\n")
                    tab_score.append(float(score))

                # On compare le score de la ligne avec le seuil
                if float(score)>score_seuil:
                    nom_fichier = line.split()[0][len(chemin):-10]
                    if nom_fichier[0]=="/":
                        nom_fichier=nom_fichier[1:]
                    chaine = nom_fichier[-3]

                    # On prend la liste de tous les résidus avec au moins un atome dans le site
                    with open(f"../data/binding_sites/{nom_fichier}.pdb", 'r') as bs:
                        for ligne in bs:
                            if ligne[21] == chaine and ligne[22:27] not in tableau_Residu:
                                tableau_Residu.append(ligne[22:27])
                                
                    Recuperation_peptide(nom_fichier, chaine, tableau_Residu, score)



# Avec les numéros de résidus, on récupère les boucles d'acides aminés consécutifs selon les paramètres
def Recuperation_peptide(nom, chaine, tab, score):
    peptide = []
    liste_peptide = []

    for i in range(0, len(tab)):

        if peptide == []:
            # Si peptide vide on ajoute le résidu considéré
            peptide.append(tab[i])

        elif int(tab[i][:-1])-int(tab[i-1][:-1]) <= gap+1:
            # Sinon si le résidu suivant a un numéro à moins de gap+1 il est ajouté
            peptide.append(tab[i])

        else:
            # Sinon si le peptide fait au moins la longueur minimal voulue alors on le conserve
            if len(peptide) >= longueur_min:
                # Pour être précis on ne récupère que les numéros du premier et du dernier résidu
                del peptide[1:-1]
                liste_peptide.append(peptide)
            peptide = []
            peptide.append(tab[i])

    # Même test de longueur mais pour le dernier peptide du tableau
    if len(peptide) >= longueur_min:
        del peptide[1:-1]
        liste_peptide.append(peptide)

    # Si au moins un peptide correpond aux paramètres, on l'écrit dans un autre fichier
    if liste_peptide != []:
        Ecriture_boucle(nom, chaine, liste_peptide, score)


# Toutes les boucles de tous les sites sont écrites dans le même fichier fasta
def Ecriture_boucle(nom, chaine, liste_peptide, score):
    with open(f"../data/dataset/{nom[:-4]}.pdb", 'r') as complexe, \
            open(f"../res/{target}_loops.fasta", 'a') as bs:
        
        for pep in liste_peptide:
            premier = pep[0]
            dernier = pep[1]
            liste_residu = ""
            liste_aa=[]
            milieu = False
            
            # On va récupérer les codes à trois lettres de la boucle avec les positions de début et de fin
            for line in complexe:
                if re.search("^(ATOM)", line) and line[21] == chaine and line[22:27] == premier and line[22:27] not in liste_aa:
                    #Premier résidu
                    liste_residu+=line[17:20]
                    liste_aa.append(line[22:27])
                    milieu = True

                elif milieu == True and line[22:27] not in liste_aa and line[22:27] != dernier:
                    # Tous les résidus au milieu
                    liste_residu+=line[17:20]
                    liste_aa.append(line[22:27])

                elif milieu == True and line[22:27] == dernier:
                    # Dernier résidu
                    liste_residu+=line[17:20]
                    break
            
            # Permet d'enlever les espaces dans les numéros de résidus
            num = ""
            for c in premier:
                if c != " ":
                    num += c
            premier = num
            num = ""
            for c in dernier:
                if c != " ":
                    num += c
            dernier = num
                    
            # Ecriture de la ligne avec le complexe, la chaine, la position et le score
            bs.write(f">{nom[:-6]}_{nom[-1]}_{premier}_{dernier}_{score}\n")
            # Ecriture de la boucle avec seq1 qui transforme les codes à trois lettres en code à une lettre
            bs.write(seq1(liste_residu)+"\n")


Recuperation_residu(target)
