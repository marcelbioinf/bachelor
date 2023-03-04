import os

files = []
for directory in os.listdir():
    if os.path.isdir(f'./{directory}'):
        if directory.endswith('idea'):
            continue
        for file in os.listdir(f'./{directory}'):
            if file.endswith('.fna'):
                files.append(f'{directory}/{file}')


file_ids = open('ZZWYNIKI/POPRAWKI/ID_LENv1.txt', 'r')
line = file_ids.readline().rstrip()

ids = {}
while line:
    seq_id, length = line.split()[0], line.split()[1]
    ids[length] = seq_id
    line = file_ids.readline().rstrip()

file_ids.close()

def get_keys(val):
    all_keys = []
    for key, value in ids.items():
        if val == value:
            all_keys.append(key)
    return all_keys

file_out = open('ZZWYNIKI/POPRAWKI/all_trnas_genomic_seqs_v1_50downup_nonrev_bothstrands.txt', 'w')

for fil in files:
    with open(fil) as file:
        for sequence in file.read().split('>'):
            if sequence:
                header, *content = sequence.splitlines()
                seqname = header[:header.find(' ')]
                if seqname in set(ids.values()):
                    all_keys = get_keys(seqname)
                    for key in all_keys:
                        positions = key.split('-')
                        if int(positions[0]) > int(positions[1]):
                            if int(positions[1]) < 50:
                                seq = ''.join(content)[(int(positions[1]) - int(positions[1])):(int(positions[0]) + 50)] #100 lub 50
                                file_out.write(
                                    f">{seqname} {int(positions[1]) - int(positions[1])}:{int(positions[0]) + 50} {(int(positions[0]) + 50) - (int(positions[1]) - int(positions[1]))}\n"
                                    f"{seq}\n")
                            else:
                                seq = ''.join(content)[(int(positions[1]) - 51):(int(positions[0]) + 50)]
                                file_out.write(
                                    f">{seqname} {int(positions[1]) - 51}:{int(positions[0]) + 50} {(int(positions[0]) + 50) - (int(positions[1]) - 51)}\n"
                                    f"{seq}\n")
                        else:
                            if int(positions[0]) < 50:
                                file_out.write(
                                    f">{seqname} {int(positions[0]) - int(positions[0])}:{int(positions[1]) + 50} {(int(positions[1]) + 50) - (int(positions[0]) - int(positions[0]))}\n"
                                    f"{''.join(content)[(int(positions[0]) - int(positions[0])):(int(positions[1]) + 50)]}\n")
                            else:
                                file_out.write(
                                    f">{seqname} {int(positions[0]) - 51}:{int(positions[1]) + 50} {(int(positions[1]) + 50) - (int(positions[0]) - 51)}\n"
                                    f"{''.join(content)[(int(positions[0]) - 51):(int(positions[1]) + 50)]}\n")
    file.close()
file_out.close()