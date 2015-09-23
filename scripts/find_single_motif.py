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
                            alignment, seq['seq'] - 1, i - 1)
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


def count_motif_occurences(alignment, motifs):
    with open(motifs) as a:
        m_list = ['motif_' + i for i in a.read().splitlines()]
    outtxt = []
    for m in m_list:
        count = 0
        for seq in alignment:
            for res in seq['sequence']:
                if m in res['features']:
                    count += 1
                    break
        outtxt.append("{} {}".format(m, count))
    return outtxt


def find_conserved_motifs(in_fasta, in_cfg, out_name, motifs):
    with open(in_fasta) as a:
        infile = unwrap(a.read().splitlines())
    with open(in_cfg) as a:
        config = a.read()
    alignment = process_alignment(infile)
    parse_config(alignment, config)
    outtxt = count_motif_occurences(alignment, motifs)
    out = open(out_name, 'w')
    out.write('\n'.join(outtxt))
    out.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("in_fasta")
    parser.add_argument("in_cfg")
    parser.add_argument("motifs")
    parser.add_argument("out")
    args = parser.parse_args()
    find_conserved_motifs(args.in_fasta, args.in_cfg, args.out,
                          args.motifs)
