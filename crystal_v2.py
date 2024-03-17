import numpy as np

class Material:
    def __init__(self, vector1, vector2, vector3, coordinates, names):
        self.A = vector1  # Basis vector a1, a2, a3 as columns of the matrix A
        self.atomic_positions = coordinates # coordinates in units of crystal basis vectors
        self.atomic_names = names

def genAtomicPositions(layer_numbers, materials):
    all_positions = []  # blank list
    total_number_layers = sum(layer_numbers)
    start_z = 0
    for material_index in materials:  # iterate thru each material
        material=materials[material_index]
        
        # I will now do a manipulation to strain this layer according to the 0th layer in the stack
        volume_material = np.det(material.A) # gives the volume of the unit cell
        strain_x = 1 - np.norm(material[0].A[:,1])/np.norm(material[material_index].A[:,1])
        strain_y = 1 - np.norm(material[0].A[:,2])/np.norm(material[material_index].A[:,2])
        print("Strain in x = " + strain_x*100 + " %, strain in y = ")
        # now, we will fix the other two 
        material[material_index].A[:,1] = material[0].A[:,1]
        material[material_index].A[:,1] = material[0].A[:,1]
        
        for layer_index in range(0, layer_numbers[material_index]-1):
            for i in range(len(material.atomic_names)):  # iterate thru each atom in each material
                name = material.atomic_names[i]
                Ri_relative = material.atomic_positions[i]  # get relative coordinates R of the ith atom
                # r = pa1 + qa2 + sa3
                Ri_absolute = material.A*Ri_relative + start_z*[0, 0, 1]
                all_positions.append({"Name": name, "Position": Ri_absolute})  # add atom to list
            start_z += A[2,2] # add on the z component of the crystal lattice basis vector to stack up the supercell in z
    return all_positions

# Example stack of Iron and Silicon Nitride
Iron = Material([2.866, 0, 0], [0, 2.866, 0], [0, 0, 2.866], [[0, 0, 0]], ['Fe'])  # 1 atom in the unit cell

# Silicon Nitride; Si-N bonds are 1.89 angstroms (materialsproject.org)
SiliconNitride = Material([1.89, 0, 0], [0, 1.89, 0], [0, 0, 1.89], [[0.25, 0.25, 0.25], [0.5, 0.5, 0.5]], ['Si', 'N'])  # example using a compound

materials = [Iron, SiliconNitride]
list_of_atoms = genAtomicPositions([1, 2], materials)  # 2 layers with two materials

# Output list
for atom in list_of_atoms:
    print(f"{atom['Name']} : {atom['Position']}")
