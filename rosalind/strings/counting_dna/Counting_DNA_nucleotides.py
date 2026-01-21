
def count_dna(word: str):
    """
    word = "AGTCGTA"
    """

    nucleo: dict = {
        'A': 0,
        'C': 0,
        'G': 0,
        'T': 0
    }

    for letter in word:
        # print(letter)
        if letter in nucleo:
            nucleo[letter] += 1
    
    return nucleo


if __name__ == '__main__':

    filename = 'rosalind_dna.txt'
    with open(filename, 'r') as file:
        contents = file.read()
        counts = count_dna(contents.strip())
        # counts.
        output = ' '.join(map(str, list(counts.values())))
        print(f'{output}')
