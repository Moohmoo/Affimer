import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys
import pymol
from Bio import Align

pymol.finish_launching()

# Read User Input
target = sys.argv[1]
query = sys.argv[2]
spath = "../test/affimer/"+query
sname = spath.split('/')[-1].split('.')[0]
chain = sys.argv[3]

# Load Structures
pymol.cmd.load(spath, sname)
pymol.cmd.disable("all")
pymol.cmd.enable(sname)

# Getting the sequence and the positions involved in the interaction with its native protein
ask_select = '(chain ' + chain + ')'
pymol.cmd.select("chain_select", ask_select)
pymol.cmd.select("bs", "not resn HOH and chain C within 5.0 of chain A")
myspace = {'positions': []}
pymol.cmd.iterate('(bs)', 'positions.append(resv)', space=myspace)
pos = []
[pos.append(position) for position in myspace['positions'] if position not in pos]

positions = []
start = 0
end = 0
for index, value in enumerate(pos) :
    if index < len(pos)-1 :
        if pos[index + 1] > value + 1 :
            end = index + 1
            positions.append(pos[start:end])
            start = end
    else :
        positions.append(pos[start:len(pos)])

sequence = pymol.cmd.get_fastastr(ask_select).split("\n")[1]
for i in range(len(positions)-1, 0, -1) :
    if (len(positions[i]) < 4) :
        positions.pop(i)

# Getting best loop
best_scores = [0 for i in range(len(positions))]
best_loops = ["" for i in range(len(positions))]
with open("../res/"+target.upper()+"_loops.fasta", "r") as filein : 
    lines = filein.readlines() 
    for i, line in enumerate(lines) : 
        if line.startswith(">") :
            score = float(line.strip().split("_")[-1])
            for j, value in enumerate(best_scores) :
                if score > value and score not in best_scores and lines[i + 1].strip() not in best_loops and len(lines[i + 1].strip()) >= 3 :
                    if abs(len(positions[j]) - 6) < len(lines[i + 1].strip()) < abs(len(positions[j]) + 6) :
                        best_scores[j] = score
                        best_loops[j] = lines[i + 1].strip()
                        break

# Change positions
myspace = {'positions': []}
pymol.cmd.iterate('(bs)', 'positions.append(resv)', space=myspace)
first_pos = myspace['positions'][0]
new_sequence = ""
for i in range(len(best_loops)) :
    pos1 = positions[i][-1] + 1 - first_pos
    if (i != len(best_loops) - 1) :
        pos2 = positions[i + 1][0] - 1 - first_pos
    else :
        pos2 = len(sequence) - 2
    
    if (pos1 >= pos2) :
        break
    new_sequence += best_loops[i]
    new_sequence += sequence[pos1 : pos2]


# Alignment global
aligner = Align.PairwiseAligner()
aligner.mode = "global"
alignments = aligner.align(sequence, new_sequence)
alignments = list(alignments[0])

# Preparing files Modeller
query_infos = f">P1;{query}\nstructureX:{query}:{first_pos}:{chain}:{len(sequence)+2}:C:DHHdom:HOMO SAPIENS : 0.0:0.0\n"
hybrid_infos = f">P1;chA\nsequence:chA: : : : :chA :: :\n"
infos = [query_infos, hybrid_infos]
with open("../res/modele.ali", "w") as fileas :
    for i in range(len(infos)) :
        fileas.write(infos[i])
        fileas.write(alignments[i])
        fileas.write("\n")
        fileas.write("*\n")
        fileas.write("\n")

# Save as png
pymol.cmd.load('pymol.pse')
pymol.cmd.png(sys.argv[1].split(".")[0]+".png")

print(f"Old affimer : {sequence}")
print(f"New affimer : {new_sequence}")

# Get out!
pymol.cmd.quit()