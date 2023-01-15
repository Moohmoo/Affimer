#!/bin/bash

export patchsearch_HOME=./../Patch_Search

MAX_POOL_SIZE=$((`nproc`/2))

function calcul_distance() {
    ./distance $1 H 5
    ./distance $1 L 5
}

function path_search() {
    chain=$1
    read_file=($2)
    target=${read_file[0]}
    query=${read_file[1]}
    antigen=${read_file[2]}
    i=${read_file[3]}

    # The script that annotates binding sites
    Rscript ../Patch_Search/SimpleAtomTyping.R ../data/binding_sites/${query}.pdb -o ../Patch_Search/Sites_Modifies/${query}_typed.pdb > /dev/null 2>&1
    # Structural comparison between the site tested and the target protein with or without a specified chain
    if [ "$chain" == "" ]; then
        Rscript ../Patch_Search/BSS.R ../Patch_Search/Proteine_Cible/${target}.pdb ../Patch_Search/Sites_Modifies/${query}_typed.pdb --bchain ${antigen} --prefix ${target} > /dev/null 2>&1
    else
        Rscript ../Patch_Search/BSS.R ../Patch_Search/Proteine_Cible/${target}.pdb ../Patch_Search/Sites_Modifies/${query}_typed.pdb --tchain ${chain} --bchain ${antigen} --prefix ${target} > /dev/null 2>&1
    fi
}

function exitfn () {
    trap SIGINT
    rm $TMP
    echo; echo -e "\n\e[1m[Commands interrupted]\e[0m"
    exit
}

trap "exitfn" INT  

export -f calcul_distance
export -f path_search

# Identification of the target protein and valid arguments
TARGET=$1
CHAIN=$2
if (( $# < 1 )) || [[ ! -f "../test/target/${TARGET}.pdb" ]]; then
    echo -e "\e[1m[ERROR]\e[0m Incorrect number of arguments or invalid arguments."
    exit 1
fi

# Directory src cleaning
make -f distance.mk clean
make -f boucles.mk clean

# Compiling the C files for the distance calculation
make -f distance.mk
start_time=$(date +%s)

# Generation of binding sites
rm -fr ../data/*/README.md
echo 'Number of cores used:' $MAX_POOL_SIZE
if [ -z "$(ls -A -- "../data/binding_sites")" ]; then
    echo "Generation of binding sites ..."
    ls -1 ../data/dataset | parallel --bar -j $MAX_POOL_SIZE calcul_distance {}
    echo -e "Binding sites data created. \xE2\x9C\x94"
else
    if [[ "$(ls -A -- "../data/binding_sites")" == *".pdb" ]]; then
        echo -e "Existing binding sites data. \xE2\x9C\x94"
    fi
fi
#ls -1 ../data/dataset | xargs -I{} bash -c 'calcul_distance "$1"'-- {}

# Execute Path Search
cp ../test/target/${TARGET}.pdb ../Patch_Search/Proteine_Cible/
TMP=tmp
i=0
if [ -z "$(ls -A -- "../res/${TARGET^^}.out")" ]; then
    echo "Executing Path Search ..."
    basename -s .pdb ../data/binding_sites/* | while read x; do echo ${TARGET^^} ${x} ${x:${#x}-3:1} $i; let "i+=1"; done > $TMP
    cat $TMP | parallel --bar -j $MAX_POOL_SIZE path_search $CHAIN {}
    rm $TMP
    mv -f *_apo.pdb ../res > /dev/null 2>&1
    mv -f *.out ../res > /dev/null 2>&1
    echo -e "Calculation of structural similarity completed. \xE2\x9C\x94"
else
    echo -e "Data on comparisons by existing Patch Search. \xE2\x9C\x94"
fi

# Getting the loops
echo "Getting the loops ..."
python3 boucles.py -t ${TARGET^^}
echo -e "Loops created. \xE2\x9C\x94"

# Preparing Modeller Files
AFFIMER=$3
AFFIMER_CHAIN=$4
python3 modeller_management.py $TARGET $AFFIMER.pdb $AFFIMER_CHAIN

# Elapsed time display
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "Elapsed time: $(date -ud "@$elapsed" +'%H hr %M min %S sec.')"

trap SIGINT
