class Material:
    def __init__(self, vector1, vector2, vector3, coordinates, names):
        self.a1 = vector1  # Basis vector a
        self.a2 = vector2  # Basis vector b
        self.a3 = vector3  # Basis vector c
        self.atomic_positions = coordinates
        self.atomic_names = names

def genAtomicPositions(layer_numbers, materials):
    all_positions = []  # blank list

    for material in materials:  # iterate thru each material
        for i in range(len(material.atomic_names)):  # iterate thru each atom in each material
            name = material.atomic_names[i]
            p, q, s = material.atomic_positions[i]  # get coordinates by index

            # r = pa1 + qa2 + sa3
            r = [
                p * material.a1[0] + q * material.a2[0] + s * material.a3[0],
                p * material.a1[1] + q * material.a2[1] + s * material.a3[1],
                p * material.a1[2] + q * material.a2[2] + s * material.a3[2]
            ]
            all_positions.append({"Name": name, "Position": r})  # add atom to list

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
