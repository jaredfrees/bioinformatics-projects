

def hamming_distance(s, t):
    count = 0
    for i in range(len(s)):
        if s[i] != t[i]:
            count += 1
    
    return count

filename = 'rosalind_hamm.txt'
with open(filename, 'r') as file:
    s = file.readline().strip()
    t = file.readline().strip()
    count = hamming_distance(s, t)
    print(s)
    print(t)

    print(count)