
seqs = {}
seq_struct = []

with open("./ZZWYNIKI/unique_trnasv2_nonrev.txt") as file:
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

out = open('t_d_loopsv3.txt', 'w')

for key in seqs.keys():
    d_start_pos = seqs[key][1].find('.......')
    if d_start_pos == -1:
        print(d_start_pos)
        continue
    d_end_pos = seqs[key][1].find('<<')
    out.write(f">{key}\nD-loop: {seqs[key][0][d_start_pos:d_end_pos]}\n")
    #print(d_start_pos, d_end_pos)
    #print(seqs[key][0][d_start_pos:d_end_pos])
    #2 rozwiazanie
    index_of_t_loop = seqs[key][1].rindex('.......')
    end_index = index_of_t_loop + 7
    if seqs[key][1][end_index-1] != '.':
        print('ERROR', seqs[key][1][end_index-1])
        continue
    out.write(f"T-loop: {seqs[key][0][index_of_t_loop:end_index]}\n")
    #print(index_of_t_loop, end_index)
    #print(seqs[key][0][index_of_t_loop:end_index])
    #1 opcja!!!
    # c_loop_toend = seqs[key][1][-20:]
    # ind = -20
    # c_end_pos = c_loop_toend.find('<')
    # if c_loop_toend[:c_end_pos] != '.......':
    # #if len(c_loop_toend[:c_end_pos]) != 7:
    #     c_loop_toend = seqs[key][1][-21:]
    #     ind = -21
    #     c_end_pos = c_loop_toend.find('<')
    #     #print(f"{key}: {c_loop_toend[:c_end_pos]}")
    #     if '.......' not in c_loop_toend[:c_end_pos]:
    #         #print(key, c_loop_toend[:c_end_pos])
    #         continue
    # out.write(f"T-loop: {seqs[key][0][ind:ind+7]}\n\n")

out.close()
