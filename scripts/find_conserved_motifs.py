#!/usr/bin/python
import argparse
import re

from jsonlibconfig import decoder


ELM_DB = "/home/joanna/data/elm_complete.txt"


def parse_elm_db():
    with open(ELM_DB) as a:
        slim_classes = a.read().splitlines()
    elm_motifs = dict()
    for i in slim_classes:
        line_i = i.split()
        elm_id = 'motif_' + line_i[0]
        elm_motifs[elm_id] = {"prob": float(line_i[2]), "GO": line_i[3:],
                              "regex": line_i[1],
                              "comp_reg": re.compile(line_i[1])}
    return elm_motifs


def make_string(sequence):
    residues = []
    for i in sequence['sequence']:
        residues.append(i['aa'])
    return ''.join(residues)


# takes 0-based seiq number and position
# position - position in the sequence without gaps
# returns position in the aligned sequence
def get_real_pos(alignment, seq, position):
    counter = 0
    real_pos = -1
    aln_sequence = alignment[seq]['sequence']
    for i in range(len(aln_sequence)):
        if aln_sequence[i]['aa'] != '-':
            counter += 1
        if position == counter - 1:
            real_pos = i
            break
    return real_pos


def predict_motifs(sequence, elm_db):
    motif_set = set()
    for elm_id in elm_db.keys():
        slim_j = elm_db[elm_id]
        reg = slim_j["comp_reg"]
        for match in reg.finditer(sequence):
            m_sp = (match.span())
            start = m_sp[0] + 1
            end = m_sp[1] + 1
            motif_set.add("{}_{}_{}".format(elm_id, start, end))
    return motif_set


def find_new_motifs(alignment, elm_db, mutations):
    sequence = re.sub('-', '', make_string(alignment[0]))
    orig_motifs = predict_motifs(sequence, elm_db)
    outtxt = []
    for mut in mutations:
        mut_seq = mutate(sequence, mut)
        mut_motifs = predict_motifs(mut_seq, elm_db)
        new = mut_motifs.difference(orig_motifs)
        if new:
            outtxt.append("new motifs on position {}".format(mut['pos']))
            outtxt.append(' '.join(list(new)))
    return outtxt


def parse_config(alignment, config):
    config = re.sub('\[', '(', re.sub('\]', ')', config))
    conf_dict = decoder.loads(config)['feature_settings']['usr_features']
    for feature in conf_dict:
        if feature['positions'] != 'empty':
            # stupid jsonlibconfig turns a list with one element into that
            # element - e.g. list with one int becomes an int,
            # so the type checks are to handle these situations
            if type(feature['positions']) == dict:
                positions = [feature['positions']]
            elif type(feature['positions']) == list:
                positions = feature['positions']
            for seq in positions:
                if seq['pos'] != 'empty':
                    if type(seq['pos']) == list:
                        pos = seq['pos']
                    elif type(seq['pos']) == int:
                        pos = [seq['pos']]
                    for i in pos:
                        real_pos = get_real_pos(
                            alignment, int(seq['seq']) - 1, int(i) - 1)
                        # add this feature to sequence number seq['seq'] -1
                        # on position real_pos
                        alignment[seq['seq'] - 1][
                            'sequence'][real_pos]['features'].append(
                            feature['name'])


def unwrap(alignment):
    new = []
    for i in alignment:
        if i.startswith('>'):
            new.append(i)
            new.append("")
        else:
            new[-1] += i
    return new


def process_alignment(infile):
    alignment = []
    seq_dict = {}
    for i in range(0, len(infile), 2):
        seq_dict['header'] = infile[i]
        seq_dict['sequence'] = [{'aa': res_j, 'features': []}
                                for res_j in infile[i + 1]]
        alignment.append(seq_dict)
        seq_dict = {}
    return alignment


def write_html():
    pass


def mutate(sequence, mutation):
    new_sequence = list(sequence)
    new_sequence[mutation['pos']] = mutation['new_aa']
    return ''.join(new_sequence)


def check_if_matches(sequence, mutation, motif_id, elm_db):
    regex = elm_db[motif_id]["comp_reg"]
    result = False
    for match in regex.finditer(sequence):
        m_sp = (match.span())
        start = m_sp[0]
        end = m_sp[1]
        if mutation['pos'] >= start and mutation['pos'] < end:
            result = True
            break
    return result


def get_motifs_on_pos(alignment, position):
    motif_ids = set()
    for i in alignment:
        for f in i['sequence'][position]['features']:
            if f.startswith('motif'):
                motif_ids.add(f)
    return list(motif_ids)


def calc_motif_conservation(alignment, mutations, elm_db):
    outtxt = []
    first_seq_string = re.sub('-', '', make_string(alignment[0]))
    for mut in mutations:
        new_seq = mutate(first_seq_string, mut)
        outtxt.append("Position {}; old aa: {}; new aa: {}".format(
            mut['pos'], first_seq_string[mut['pos']], mut['new_aa']))
        outtxt.append("motif id; cons_all; cons_100; cons_50; original seq; mut seq")
        real_pos = get_real_pos(alignment, 0, mut['pos'])
        # motif counts holds counts for all occurences of the motif,
        # occurences in the 1st 50 sequnces, and in the first 100 sequences
        motif_ids = get_motifs_on_pos(alignment, real_pos)
        # first_motifs = [i for i in alignment[0]['sequence'][real_pos]['features'] if i.startswith('motif')]
        # motif_counts = {feat: {'all': 1, '50': 1, '100': 1}
        #                 for feat in first_motifs}
        motif_counts = {feat: {'all': 0, '50': 0, '100': 0}
                        for feat in motif_ids}

        for i, seq in enumerate(alignment):
            for feat in motif_counts:
                found = False
                for pos in range(real_pos - 10, real_pos + 11):
                    if (pos > -1 and pos < len(seq['sequence'])
                            and feat in seq['sequence'][pos]['features']):
                        found = True
                        break
                # for feat in set(res['features']):
                #     if feat not in motif_counts.keys():
                #         motif_counts[feat] = {'all': 0, '50': 0, '100': 0}
                if found:
                    motif_counts[feat]['all'] += 1
                    # i < 49 because first sequence is already included
                    if i < 49:
                        motif_counts[feat]['50'] += 1
                    if i < 99:
                        motif_counts[feat]['100'] += 1
        for m in motif_counts:
            n_50 = 50
            n_100 = 100
            if len(alignment) < 50:
                n_50 = len(alignment)
            if len(alignment) < 100:
                n_100 = len(alignment)
            # motif_counts[m] = float(motif_counts[m]) / len(alignment)
            c_all = round(float(motif_counts[m]['all']) / len(alignment), 2)
            c_50 = round(float(motif_counts[m]['50']) / n_50, 2)
            c_100 = round(float(motif_counts[m]['100']) / n_100, 2)
            if m.startswith('motif_'):
                present_in_old = check_if_matches(first_seq_string,
                                                  mut, m, elm_db)
                present_in_new = check_if_matches(new_seq, mut, m, elm_db)
            else:
                present_in_old = ""
                present_in_new = ""
            # if any([c_all >= 0.05, c_100 >= 0.05,
            #         c_50 >= 0.05, present_in_old]):
            # TEMPORARY
            # to output all features appearing on this position
            if present_in_old or present_in_new:
                outtxt.append(
                    "{}; {}; {}; {}; {}; {}".format(m, c_all, c_100, c_50,
                                                    present_in_old,
                                                    present_in_new))
    return outtxt


def calc_ptm_conservation(alignment, mutations):
    outtxt = []
    first_seq_string = re.sub('-', '', make_string(alignment[0]))
    for mut in mutations:
        ptm_counts = {}
        outtxt.append("Position {}; old aa: {}; new aa: {}".format(
            mut['pos'], first_seq_string[mut['pos']], mut['new_aa']))
        outtxt.append("ptm type; cons_all; cons_100; cons_50; original seq;")
        for i, seq in enumerate(alignment):
            real_pos = get_real_pos(alignment, 0, mut['pos'])
            res = seq['sequence'][real_pos]
            feat_set = set()
            for feat in set(res['features']):
                if feat.startswith('ptm'):
                    # count all annotation levels as 1
                    if feat != "ptm_phosphP":
                        ptm_key = feat[:-1]
                    else:
                        ptm_key = feat
                    if ptm_key not in feat_set:
                        feat_set.add(ptm_key)
                        if ptm_key not in ptm_counts.keys():
                            ptm_counts[ptm_key] = {'all': 0, '50': 0, '100': 0}
                        ptm_counts[ptm_key]['all'] += 1
                        # i < 49 because first sequence is already included
                        if i < 50:
                            ptm_counts[ptm_key]['50'] += 1
                        if i < 100:
                            ptm_counts[ptm_key]['100'] += 1
        for m in ptm_counts:
            n_50 = 50
            n_100 = 100
            if len(alignment) < 50:
                n_50 = len(alignment)
            if len(alignment) < 100:
                n_100 = len(alignment)
            # motif_counts[m] = float(motif_counts[m]) / len(alignment)
            c_all = round(float(ptm_counts[m]['all']) / len(alignment), 2)
            c_50 = round(float(ptm_counts[m]['50']) / n_50, 2)
            c_100 = round(float(ptm_counts[m]['100']) / n_100, 2)
            check_first = [i.startswith(ptm_key)
                           for i in alignment[0]['sequence'][real_pos]['features']]
            present_in_old = any(check_first)
            outtxt.append(
                "{}; {}; {}; {}; {};".format(m, c_all, c_100, c_50,
                                             present_in_old))
    return outtxt


'''


def calc_motif_conservation(alignment, mutations, elm_db):
    outtxt = []
    first_seq_string = re.sub('-', '', make_string(alignment[0]))
    for mut in mutations:
        new_seq = mutate(first_seq_string, mut)
        outtxt.append("Position {}; old aa: {}; new aa: {}".format(
            mut['pos'], first_seq_string[mut['pos']], mut['new_aa']))
        outtxt.append("motif id; cons_all; cons_100; cons_50; original seq; mut seq")
        real_pos = get_real_pos(alignment, 0, mut['pos'])
        # motif counts holds counts for all occurences of the motif,
        # occurences in the 1st 50 sequnces, and in the first 100 sequences

        # motif_counts = {feat: {'all': 0, '50': 0, '100': 0}for feat in
        #                 alignment[0]['sequence'][real_pos]['features']}
        motif_counts = {}
        for i, seq in enumerate(alignment):
            res = seq['sequence'][real_pos]
            for feat in set(res['features']):
                if feat not in motif_counts.keys():
                    motif_counts[feat] = {'all': 0, '50': 0, '100': 0}
                motif_counts[feat]['all'] += 1
                # i < 49 because first sequence is already included
                if i < 50:
                    motif_counts[feat]['50'] += 1
                if i < 100:
                    motif_counts[feat]['100'] += 1
        for m in motif_counts:
            # motif_counts[m] = float(motif_counts[m]) / len(alignment)
            n_50 = 50
            n_100 = 100
            if len(alignment) < 50:
                n_50 = len(alignment)
            if len(alignment) < 100:
                n_100 = len(alignment)
            c_all = round(float(motif_counts[m]['all']) / len(alignment), 2)
            c_50 = round(float(motif_counts[m]['50']) / n_50, 2)
            c_100 = round(float(motif_counts[m]['100']) / n_100, 2)
            if m.startswith('motif_'):
                present_in_old = check_if_matches(first_seq_string,
                                                  mut, m, elm_db)
                present_in_new = check_if_matches(new_seq, mut, m, elm_db)
            else:
                present_in_old = ""
                present_in_new = ""
            if any([c_all >= 0.05, c_100 >= 0.05,
                    c_50 >= 0.05, present_in_old, m.startswith('ptm')]):
                outtxt.append(
                    "{}; {}; {}; {}; {}; {}".format(m, c_all, c_100, c_50,
                                                    present_in_old,
                                                    present_in_new))
    return outtxt
'''


def find_conserved_motifs(in_fasta, in_cfg, out_name, in_mutations):
    with open(in_fasta) as a:
        infile = unwrap(a.read().splitlines())
    with open(in_cfg) as a:
        config = a.read()
    with open(in_mutations) as a:
        mutations = [{'pos': int(i.split()[0]), 'new_aa': i.split()[1]}
                      for i in a.read().splitlines()]
    elm_db = parse_elm_db()
    alignment = process_alignment(infile)
    parse_config(alignment, config)
    outtxt = calc_motif_conservation(alignment, mutations, elm_db)
    outtxt.extend(calc_ptm_conservation(alignment, mutations))
    outtxt.extend(find_new_motifs(alignment, elm_db, mutations))
    out = open(out_name, 'w')
    out.write('\n'.join(outtxt))
    out.close()
    # write_html(alignment_dict, out_html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("in_fasta")
    parser.add_argument("in_cfg")
    parser.add_argument("in_mutations")
    parser.add_argument("out")
    args = parser.parse_args()
    find_conserved_motifs(args.in_fasta, args.in_cfg, args.out,
                          args.in_mutations)
