seqs = {}
seq_struct = []

with open("./ZZWYNIKI/unique_trnasv2.txt") as file:
    for line in file:
        if line.startswith('>'):
            header = line.split()[1]
        elif line.upper().startswith('SEQ: '):
            seq = line.upper().lstrip('SEQ: ').rstrip("\n")
            seq_struct.append(seq)
        elif line.upper().startswith('STR: '):
            struct = line.upper().lstrip('STR: ').rstrip('\n')
            seq_struct.append(struct)
        elif line == "\n":
            seqs[header] = seq_struct.copy()
            seq_struct.clear()
file.close()

out = open('ZZWYNIKI/box_and_loops.txt', 'w')
ignore = []
for key in seqs.keys():
    a_start_pos = seqs[key][1].find('..')
    if a_start_pos == 7 or a_start_pos == 6:
        d_end_pos = seqs[key][1].find('<<')
        abox_dloop = seqs[key][0][a_start_pos:d_end_pos]
    else:
        ignore.append(seqs[key][0])
        continue
    out.write(f">{key}\nAboxDlopp: {abox_dloop}\n")
    index_of_t_loop = seqs[key][1].rindex('.......')
    index_b_box = index_of_t_loop - 1
    bbox_end_index = index_of_t_loop + 7 + 2
    out.write(f'BboxTloop: {seqs[key][0][index_b_box:bbox_end_index]}\n')

