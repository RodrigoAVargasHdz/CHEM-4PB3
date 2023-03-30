from rdkit.Chem import Lipinski, MolFromSmiles
from rdkit import Chem
from rdkit.Chem import AllChem
from evomol import run_model
import subprocess 
import numpy as np

def energy_from_lammps(smiles):

    # get the atomic positions using rdkit 
    m = Chem.MolFromSmiles(smiles)                # rdkit structure from smiles
    molecule = Chem.AddHs(m)                      # adding hydrogen atoms
    AllChem.EmbedMolecule(molecule)               # I don't know what this is doing
    AllChem.UFFOptimizeMolecule(molecule)         # first optimization with UFF 
    molecule.GetConformer()                       # I don't know what this is doing
    atomic_positions = []
    for i, atom in enumerate(molecule.GetAtoms()):
        atom_position = []
        positions = molecule.GetConformer().GetAtomPosition(i)
        atom_position.append(atom.GetSymbol())
        atom_position.append(positions.x)
        atom_position.append(positions.y)
        atom_position.append(positions.z)
        atomic_positions.append(atom_position)
    
    # change atomic symbol to letter to run lammps
    for i in range(0,len(atomic_positions)):
        if (atomic_positions[i][0]=="C"):
           atomic_positions[i][0]="1"
        elif(atomic_positions[i][0]=="H"):
            atomic_positions[i][0]="2"
        elif(atomic_positions[i][0]=="O"):
            atomic_positions[i][0]="3"
        else:
            atomic_positions[i][0]="4"

    # write atomic positions in system.data
        
        # getting the initial data from system.data

    data_line = []
    with open("system.data","r") as file:
        for line in file:
            data_line_aux = []
            line = line.split()
            for item in line:
                if (item=="X"):
                    item = len(atomic_positions) # changing X in system.data to add the number of atoms in the system
                    item += 2                    # adding two because of N2
                data_line_aux.append(item)
            data_line.append(data_line_aux)

        # adding the atoms in the possible molecule
    number_atom_lammps = 3 # it starts with 3 because in the script there already is a N2 molecule
    for atom in atomic_positions:
        data_line_aux = []
        data_line_aux.append(str(number_atom_lammps))
        number_atom_lammps += 1
        data_line_aux.append(str("2")) # number of the molecule
        data_line_aux.append(atom[0])
        data_line_aux.append("0.0000")
        data_line_aux.append(atom[1])
        data_line_aux.append(atom[2])
        data_line_aux.append(atom[3])
        data_line.append(data_line_aux)

        # writing the data to run.data 

    with open("run.data","w") as file:
        for line in data_line:
            for item in line:
                file.write(str(item))
                file.write(" ")
            file.write("\n")

    # run lammps
    subprocess.run("lmp < system_minimization.in", shell=True)

    # get affinity energy from output.out

    with open("log.lammps","r") as file:
        data = []
        data_tentative_molecule = []
        both_molecules = True # getting the total energy of the system 
        for line in file:
            line = line.split()
            try: 
                if (line[0]=="delete_atoms"):
                    both_molecules = False
            except:
                continue          
            if (len(line)==8):   
                if (both_molecules):
                    data.append(line[3])
                    data = data[-1000:]
                    # removing the words from the log.lammps file
                    for i in range(0,len(data)): 
                        try:
                            data[i] = float(data[i])
                        except:
                            data.pop(i)                       
                else:
                    data_tentative_molecule.append(line[3])
                    data_tentative_molecule = data_tentative_molecule[-1000:]
                    for i in range(0,len(data_tentative_molecule)):
                        try:
                            data_tentative_molecule[i] = float(data_tentative_molecule[i])
                        except:
                            data_tentative_molecule.pop(i)

    print("Adsorption Energy ", np.mean(data)-np.mean(data_tentative_molecule)-(-192.085635))
    return np.mean(data)-np.mean(data_tentative_molecule)-(-192.085635)



objective_tree = {
    # The objective function is a linear combination of three sub-functions
    "type": "linear_combination",
    "coef": [1],
    "functions": [
        {
            # Centering a Gaussian function on the adsorption energy of -10 kcal/mol 
            "type": "opposite",
            "function":  (energy_from_lammps, "energy_from_lammps")
        },
    ]
}


run_model({
    "obj_function": objective_tree,
    "action_space_parameters": {
        "atoms" : "C,H,O,N"
    },
    "io_parameters": {
        "model_path": "teste_rodrigo"
    },
    "optimization_parameters": {
     "max_steps": 5
    }
})
