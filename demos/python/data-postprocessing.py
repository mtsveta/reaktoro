data_folder = 'demos/python/data-files/'
file_with_species = 'species.txt'
file_with_species_codes = 'species-codes.txt'

def get_data_from_file(filepath):
    items = []
    with open(filepath,'r') as file:
        for line in file:
            for word in line.split():
                items.append(word.replace("'", "")) #items.append(re.sub(r'[^\w]', '', word))
    return items

species = get_data_from_file(data_folder + file_with_species)
species_codes = get_data_from_file(data_folder + file_with_species_codes)

print(species)
print(species_codes)

# Fetch all the aqueous species
aqueous_species = [] # codes 'S', 'T', 'W'
gaseous_species = [] # codes 'G', 'M', 'I'
minerals = [] # code 'O'
for item in zip(species, species_codes):
    if item[1] == 'S' or item[1] == 'T' or item[1] == 'W': aqueous_species.append(item[0])
    if item[1] == 'G' or item[1] == 'M' or item[1] == 'I': gaseous_species.append(item[0])
    if item[1] == 'O': minerals.append(item[0])
print(aqueous_species)
print(gaseous_species)
print(minerals)
