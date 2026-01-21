
filename = 'rna_to_dna.txt'
with open(filename, 'r') as file:
    contents = file.read().strip()
    
    output = ''
    for letter in contents:
        if letter == 'T':
            letter = 'U'
        output += letter
    
    print(output)
