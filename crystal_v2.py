import numpy as np

class Material:
    def __init__(self, vector1, vector2, vector3, coordinates, names):
        self.A = np.array([vector1,vector2,vector3]).T  # Basis vector a1, a2, a3 as columns of the matrix A
        self.atomic_positions = coordinates # coordinates in units of crystal basis vectors
        self.atomic_names = names

def genAtomicPositions(layer_numbers, materials):
    all_positions = []  # blank list
    total_number_layers = sum(layer_numbers)
    start_z = 0
    for material_index in range(len(materials)):  # iterate thru each material
        material = materials[material_index]
        
        # I will now do a manipulation to strain this layer according to the 0th layer in the stack
        volume_material = np.linalg.det(material.A)  # gives the volume of the unit cell
        strain_x = 1 - np.linalg.norm(materials[0].A[:,1])/np.linalg.norm(material.A[:,1])
        strain_y = 1 - np.linalg.norm(materials[0].A[:,2])/np.linalg.norm(material.A[:,2])
        print(f"Strain in x = {strain_x*100} %, strain in y = {strain_y*100} %")
        
        # now, we will fix the other two 
        material.A[:,1] = materials[0].A[:,1]
        material.A[:,2] = materials[0].A[:,2]
        
        for layer_index in range(layer_numbers[material_index]):
            for i in range(len(material.atomic_names)):  # iterate thru each atom in each material
                name = material.atomic_names[i]
                Ri_relative = np.array(material.atomic_positions[i])  # get relative coordinates R of the ith atom
                
                # r = pa1 + qa2 + sa3, matrix multiplication for coordinate transformation
                Ri_absolute = material.A @ Ri_relative + start_z * np.array([0, 0, 1])
                all_positions.append({"Name": name, "Position": Ri_absolute})  # add atom to list

            start_z += material.A[2,2]  # add on the z component of the crystal lattice basis vector to stack up the supercell in z

    return all_positions

#function that counts the total number of atoms in the stack as well as how many of each element
#takes the 'stack' as the input
def count_atoms(atoms_list):
    element_counts = {} #define an empty list
    for atom in atoms_list: 
        element = atom['Name']
        if element in element_counts:#check for repeats
            element_counts[element] += 1
        else:
            element_counts[element] = 1 #first occurrence, add to list

    total_atoms = sum(element_counts.values())
    print(f"{total_atoms} total atoms")
    for element, count in element_counts.items(): # print each element + its corresponding count
        print(f"{element}: {count} atom(s)")


# Example stack of Iron and Silicon Nitride
Iron = Material([2.866, 0, 0], [0, 2.866, 0], [0, 0, 2.866], [[0, 0, 0]], ['Fe'])  # 1 atom in the unit cell

# Silicon Nitride; Si-N bonds are 1.89 angstroms (materialsproject.org)
SiliconNitride = Material([1.89, 0, 0], [0, 1.89, 0], [0, 0, 1.89], [[0.25, 0.25, 0.25], [0.5, 0.5, 0.5]], ['Si', 'N'])  # example using a compound

materials = [Iron, SiliconNitride, Iron]
list_of_atoms = genAtomicPositions([5,5,5], materials)  # 2 layers with two materials

#atoms count
count_atoms(list_of_atoms)

# Output list
for atom in list_of_atoms:
    print(f"{atom['Name']} : {atom['Position']}")
    
    
    
