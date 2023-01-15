# Load standard Modeller classes.
from modeller import *
# Load the automodel class.
from modeller.automodel import *
import sys

modele = sys.argv[1]
nombre = sys.argv[2]


# Request verbose output.
log.verbose()
# Create a new MODELLER environment to build this model in.
env = environ()

# Read in HETATM records from template PDBs.
env.io.hetatm = True

a = automodel(
    env,
    alnfile="../res/modele.ali",   # Alignment filename.
    knowns=[modele],       # Codes of the templates.
    sequence="chA",               # Code of the target.
    assess_methods=(assess.DOPE, assess.GA341)
)

# Index of the first model.
a.starting_model = 1
# Index of the last model.
a.ending_model = nombre
# Do the actual homology modeling.
a.make()
