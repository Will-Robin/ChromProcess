from pathlib import Path

path = Path(__file__)
script_dir = path.parent

def readTwoColInfo(file, col_idx):
    assignments = {}
    with open(file, 'r') as f:
        for c, line in enumerate(f):
            if c > 0:
                spl = line.strip('\n').split(',')
                assignments[spl[0]] = spl[col_idx]
            else:
                pass
    return assignments

info_container = []
with open(script_dir/'compound_properties.csv', 'r') as f:
    for c, line in enumerate(f):
        if c == 0:
            header = line.strip('\n').split(',')
        else:
            spl = line.strip('\n').split(',')
            info_container.append(spl)

info_container = [list(i) for i in zip(*info_container)]

props_dict = {}
for n,i in zip(header, info_container):
    props_dict[n] = i

colour_assignments = {k:v for k,v in zip(props_dict['@ SMILES'], props_dict['colour'])}
mms = {k:float(v) for k,v in zip(props_dict['@ SMILES'], props_dict['Mr_gmol-1'])}
canonical_SMILES = {k:v for k,v in zip(props_dict['compound_name'], props_dict['@ SMILES'])}
smiles_to_names = {}
init_can_SMILES = [c for c in canonical_SMILES]
for c in init_can_SMILES:
    spl_name = c.split(' ')
    canonical_SMILES[spl_name[0]] = canonical_SMILES[c]
    smiles_to_names[canonical_SMILES[c]] = spl_name[0]

class_assignments =  {k:v for k,v in zip(props_dict['@ SMILES'], props_dict['Class'])}

fr_assign = readTwoColInfo(script_dir/'frag_assignments.csv', 1)
frag_assignments = {float(f):fr_assign[f] for f in fr_assign}
frag_colours = readTwoColInfo(script_dir/'frag_assignments.csv', 2)
