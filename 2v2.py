import os

files = []                                                                                     #BLEDY
for directory in os.listdir():
    if os.path.isdir(f'./{directory}'):
        if directory.endswith('idea'):
            continue
        for file in os.listdir(f'./{directory}'):
            if file.endswith('rss_hc.txt'):
                files.append(f'{directory}/{file}')
                                                                                            #files.append([f'{directory}/{file}' for file in os.listdir(f'./{directory}') if file.endswith('rss_hc.txt')]) - to mi daje tylko ostatni plik onelinera ine moge tak dac w petli

seq, helper, all_sequences = [], [], []
line_counter = 0                                                                            #enumerate nie można edytować pózniej
for file in files:
    seqs = open(f'./{file}', 'r')
    for line in seqs:
        if line == '\n':
            line_counter = 0
            all_sequences.append(seq.copy())
            helper.clear()
            seq.clear()
            continue
        if line_counter < 3:
            helper.extend(line.split())
            if line_counter == 2:
                if line.startswith('Possible '):
                    continue
                seq.append(tuple(helper))
        else:
            seq.append(line.rstrip('\n').upper())
        line_counter += 1

#after this I have a list of lists. Each list represents one sequence and has tuple as header and few strings as the rest

# for sequence in all_sequences:
#     pos = sequence[0][1].lstrip('(').rstrip(')').split('-')
#     if int(pos[0]) > int(pos[1]):
#         new_seq = sequence[2].lstrip('Seq: ')
#         sequence[2] = f"Seq: {new_seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))}"


all_sequences_dict = {}
for sequence in all_sequences:
    if sequence[2] in all_sequences_dict.keys():      #this let us have only best scored sequences - niepotrzebne bo tak samo punktowane xddd
        if all_sequences_dict[sequence[2]][0][13] < sequence[0][13]:
            all_sequences_dict[sequence[2]] = [item for i, item in enumerate(sequence) if i != 2]
    else:
        all_sequences_dict[sequence[2]] = [item for i, item in enumerate(sequence) if i != 2]

print(len(all_sequences_dict))   # ostateczna ilosc unikalnych tRNA  - 2546 # 2981

file_out = open('ZZWYNIKI/unique_trnasv3_nonrev.txt', 'w')

for seq in all_sequences_dict.keys():
    file_out.write((f">{all_sequences_dict[seq][0][0].split('t')[0].rstrip('.')}> {' '.join(all_sequences_dict[seq][0])}\n{all_sequences_dict[seq][1]}\n{seq}\n{all_sequences_dict[seq][2]}\n\n"))

file_out = open('ID_LENv3.txt', 'w')
for seq in all_sequences_dict.keys():
   file_out.write(f"{all_sequences_dict[seq][0][0].split('t')[0].rstrip('.')} {all_sequences_dict[seq][0][1].rstrip(')').lstrip('(')}\n")