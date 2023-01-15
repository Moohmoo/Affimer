

BSBANK=/home/concorde/bochinski/Documents/Projet/Patch_Search/BS_Modifie

export patchsearch_HOME=/home/concorde/bochinski/Documents/Projet/Patch_Search

Rscript  $patchsearch_HOME/BSS.R $2.pdb $BSBANK -p 1 -o $2.out --prefix $2

python3 boucle.py -t $2
