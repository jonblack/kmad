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
    return alignment


def write_html():
    pass

def filter_conserved_motifs(alignment):
    for res_i in


def convert_to_html(in_fasta, in_cfg, out_html, mode):
    with open(in_fasta) as a:
        infile = unwrap(a.read().splitlines())
    with open(in_cfg) as a:
        config = a.read()
    elm_db = parse_elm_db()
    alignment = process_alignment(infile)
    parse_config(alignment, config)
    filter_conserved_motifs(alignment, elm_db)
    # write_html(alignment_dict, out_html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("in_fasta")
    parser.add_argument("in_cfg")
    parser.add_argument("out_html")
    parser.add_argument("mode")
    args = parser.parse_args()
    convert_to_html(args.in_fasta, args.in_cfg, args.out_html, args.mode)
