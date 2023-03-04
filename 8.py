seqs = {}
seq_struct = []

with open("ZZWYNIKI/unique_trnasv2.txt") as file:
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

#out = open('ZZWYNIKI/a_b_box.txt', 'w')
#out = open('ZZWYNIKI/a_boxv-1_11.txt', 'w')
#out11 = open('ZZWYNIKI/a_boxv-1_12.txt', 'w')
out2 = open('ZZWYNIKI/b_boxv-1.txt', 'w')

ignore = []
for key in seqs.keys():
    a_start_pos = seqs[key][1].find('..>')
    if a_start_pos == 7:
        if seqs[key][0][a_start_pos+1] == 'T':
            a_start_pos += 1
        abox = seqs[key][0][a_start_pos:a_start_pos+11]
        if not abox.endswith('GG'):
            abox = seqs[key][0][a_start_pos:a_start_pos + 12]
    elif a_start_pos == 6:
        if seqs[key][0][a_start_pos-1] == 'T':
            a_start_pos -= 1
        abox = seqs[key][0][a_start_pos:a_start_pos+12]
        if not abox.endswith('GG'):
            abox = seqs[key][0][a_start_pos:a_start_pos + 13]
    else:
        ignore.append(seqs[key][0])
        continue
    if len(abox) > 12:
        print(abox)
    #if len(abox) == 12:
        #out11.write(f">{key}\n{abox}\n")
    #else:
        #out.write(f">{key}\n{abox}\n")
    #BBOX
    index_of_t_loop = seqs[key][1].rindex('.......')
    index_b_box = index_of_t_loop - 2  #lub - 2
    bbox_end_index = index_of_t_loop + 7 + 2 #lub +3
    out2.write(f'>{key}\n{seqs[key][0][index_b_box:bbox_end_index]}\n')

print(len(ignore))
out2.close()
