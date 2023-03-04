import os
from collections import Counter
final = [0, 0, 0, 0]

files = []
for directory in os.listdir():
    if os.path.isdir(f'./{directory}'):
        if directory.endswith('idea'):
            continue
        for file in os.listdir(f'./{directory}'):
            if file.endswith('.fna'):
                files.append(f'{directory}/{file}')
print(files)
#for fil in files:
with open(files[0]) as file:
    global_num = [0, 0, 0, 0]
    length = 0
    for chromosome in file.read().split('>'):
        if chromosome:
            header, *content = chromosome.splitlines()
            seq = ''.join(content).upper()
            numbers = dict(zip(['A', 'T', 'G', 'C'], list(map(seq.count, "ATGC"))))
            print(numbers)
            for i, value in enumerate(numbers.values()):
                global_num[i] += value
            # print(*(seq.count(nuc) for nuc in "ATGC"))  # 2 opcja
            length += sum(numbers.values())
    print(f"{header}: {[i / length for i in global_num]}\n")
    #final += [i / length for i in global_num]  -> Zła komenda
    #final = list(map(lambda x: x + x / length, global_num)) tez zla komenda

print([i / 11 for i in final])




