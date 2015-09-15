#!/usr/bin/python
import argparse
import logging
import math
import os
import re
import subprocess
import tempfile
import urllib2

from jsonlibconfig import encoder
from httplib import BadStatusLine
from socket import error as SocketError


SWISS_FASTA_DIR = "/home/joanna/data/swiss_fasta/uniprot_fasta"
SWISS_DAT_DIR = "/home/joanna/data/swissprot_dat/uniprot_dat"
SWISS_BLAST = "/home/joanna/data/uniprot_sprot"
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


def read_fasta(fastafilename):
    with open(fastafilename) as a:
        fastafile = a.read()
    fastafile_list = fastafile.splitlines()
    result = []
    for i, elementI in enumerate(fastafile_list):
        if elementI.startswith('>') or fastafile_list[i-1].startswith('>'):
            result += [elementI]
        else:
            result[-1] += elementI
    return result


def get_id(sequence):
    if len(sequence.split('|')) >= 3:
        return sequence.split('|')[2].split(' ')[0]
    else:
        return sequence.lstrip('>')


def check_id(uniprot_id, seq):
    result = False
    path = os.path.join(SWISS_FASTA_DIR, uniprot_id + '.fasta')
    uni_seq = ''
    if os.path.exists(path):
        with open(path) as a:
            uni_seq = ''.join(a.read().splitlines()[1:])
    if uni_seq == seq:
        result = True
    return result


def run_blast(sequence, blastdb):
    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    with tmp_file as f:
        f.write(sequence)
    out_blast = tmp_file.name + '.blastp'
    args = ["blastp", "-query", tmp_file.name, "-evalue", "1e-5",
            "-num_threads", "15", "-db", blastdb,
            "-out", out_blast, '-outfmt', '10',
            "-max_target_seqs", '10']
    try:
        subprocess.call(args)
    except subprocess.CalledProcessError as e:
        print "Error: {}".format(e.output)

    if os.path.exists(out_blast):
        with open(out_blast) as a:
            output = a.read().splitlines()
        logging.debug(output)
        if output:
            logging.debug('output True')
        os.remove(out_blast)
    else:
        output = []
    os.remove(tmp_file.name)
    return output


def find_seq_id_from_blast(sequence):
    blast_result = run_blast(sequence, SWISS_BLAST)
    # if blast_result and float(blast_result[0].split(',')[2]) >= 70:
    seq_id = ""
    if blast_result:
        firstline = blast_result[0].split(',')
        if (firstline[2] == "100.00"
                and firstline[4] == "0" and firstline[5] == "0"
                and firstline[10] == "0.0"
                and firstline[6:8] == firstline[8:10]):
            seq_id = firstline[1].split('|')[1]
    return seq_id


def find_uniprot_ids(fasta_seqs):
    seq_ids = []
    for i in range(0, len(fasta_seqs), 2):
        header_id = get_id(fasta_seqs[i])
        if check_id(header_id, fasta_seqs[i + 1]):
            seq_ids.append(header_id)
        else:
            seq_id = find_seq_id_from_blast(fasta_seqs[i + 1])
            seq_ids.append(seq_id)
    return seq_ids


def remove_gaps(fasta_seqs):
    fasta_seqs_degapped = []
    for i in fasta_seqs:
        # check if it's a header line
        if i.startswith('>'):
            fasta_seqs_degapped.append(i)
        else:
            fasta_seqs_degapped.append(re.sub('-', '', i))
    return fasta_seqs_degapped


def format_positions(positions):
    result = []
    for i in positions:
        result.append({"seq": i,
                       "pos": tuple(positions[i])})
    return tuple(result)


def make_conf_dict(metadata):
    conf_dict = {"feature_settings": {"usr_features": []}}
    for feature in metadata:
        if feature.startswith('motif_'):
            al_score = metadata[feature]['prob']
        else:
            al_score = 0
        positions = format_positions(metadata[feature]['positions'])
        single_feat = {"name": feature,
                       "add_score": al_score,
                       "subtract_score": 0,
                       "add_features": [],
                       "add_tags": [],
                       "add_exceptions": [],
                       "subtract_features": [],
                       "subtract_tags": [],
                       "subtract_exceptions": [],
                       "subtract_features": [],
                       "pattern": '',
                       "positions": positions}
        conf_dict["feature_settings"]["usr_features"].append(
            single_feat)
    return conf_dict


def write_conf_file(metadata, outname):
    data_dict = make_conf_dict(metadata)
    indent = 2
    outtxt = encoder.dumps(data_dict, indent=indent)
    out = open(outname, 'w')
    out.write(outtxt)
    out.close()


def get_uniprot_data(uni_id):
    data_path = os.path.join(SWISS_DAT_DIR, uni_id + '.dat')
    seq_data = ""
    if os.path.exists(data_path):
        with open(data_path) as a:
            seq_data = a.read().splitlines()
    else:
        try:
            url = "http://www.uniprot.org/uniprot/{}.txt".format(uni_id)
            req = urllib2.Request(url)
            seq_data = urllib2.urlopen(req).read().splitlines()
        except (urllib2.HTTPError,  urllib2.URLError,
                SocketError, BadStatusLine):
            print "Couldn't find a uniprot entry for id {}".format(
                uni_id)
    return seq_data


# merge dicts with structure:
# {
#   'motif_id': {
#            'prob': 1.0,
#            'positions': {1: [1, 2, 3]}
#          }
# }
#
def merge_motif_dicts(x, y):
    z = {}
    for i in set(x.keys() + y.keys()):
        z[i] = {}
        if i in x and i in y:
            z[i]['positions'] = merge_list_dicts(x[i]['positions'],
                                                 y[i]['positions'])
            z[i]['prob'] = x[i]['prob']
        elif i in x:
            z[i] = x[i]
        elif i in y:
            z[i] = y[i]
    return z


# merge dicts with structure:
# {'key': list()
# }
# keep only unique list elements
def merge_list_dicts(x, y):
    z = {}
    for i in set(x.keys() + y.keys()):
        z[i] = []
        if i in x:
            z[i].extend(x[i])
        if i in y:
            z[i].extend(y[i])
        z[i] = list(set(z[i]))
    return z


def run_netphos(filename):
    phosphorylations = set([])
    try:
        args = ['netphos', filename]
        netphos_result = subprocess.check_output(args).splitlines()

        for lineI in netphos_result:
            if len(lineI.split()) > 0 and lineI.split()[-1] == 'YES':
                phosphorylations.add(int(lineI.split()[2]))
    except subprocess.CalledProcessError:
        print "Not running netphos"
    return list(phosphorylations)


def get_ptm_type(txt_line):
    ptm_found = True
    # check ptm_type kind and insert site in the results list
    if "Phospho" in txt_line and "MOD_RES" in txt_line:
        ptm_type = 'phosph'
    elif "amide" in txt_line and "MOD_RES":
        ptm_type = 'amid'
    elif "CARBOHYD" in txt_line and "O-linked" in txt_line:
        ptm_type = 'Oglyc'
    elif "CARBOHYD" in txt_line and "N-linked" in txt_line:
        ptm_type = 'Nglyc'
    elif "MOD_RES" in txt_line and "acetyl" in txt_line:
        ptm_type = 'acetyl'
    elif "MOD_RES" in txt_line and "hydroxy" in txt_line:
        ptm_type = 'hydrox'
    elif "MOD_RES" in txt_line and "methyl" in txt_line:
        ptm_type = 'meth'
    else:
        ptm_found = False
        ptm_type = ""
    return ptm_type, ptm_found


def get_annotation_level(uni_features):
    levels_dict = {'269': 0, '314': 0, '353': 0, '315': 0, '316': 0, '270': 0,
                   '250': 1, '266': 1, '247': 1, '255': 1, '317': 1, '318': 1,
                   '319': 1, '320': 1, '321': 1, '245': 1,
                   '304': 2, '303': 2, '305': 2, '307': 2,
                   '501': 3}
    i = 0
    n = 3
    reading = True
    while reading and i < len(uni_features):
        if uni_features[i].startswith("FT          ") or i == 0:
            if "ECO:0000" in uni_features[i]:
                start = uni_features[i].index("ECO:0000") + 8
                eco_code = uni_features[i][start:start+3]
                n = levels_dict[eco_code]
                break
            i += 1
        else:
            reading = False
    return n


# get PTM annotations from Uniprot's sequence data
# 0-based!!
def parse_seq_data(seq_data, ptm_dict, seq_index):
    for i, lineI in enumerate(seq_data):
        if (lineI.startswith('FT')
                and len(lineI.split()) > 3
                and lineI.split()[3].isdigit()):
            ptm_type, ptm_found = get_ptm_type(lineI)
            if ptm_found:
                n = str(get_annotation_level(seq_data[i:i + 10]))
                position = int(lineI.split()[3])
                ptm_name = 'ptm_' + ptm_type + n
                if seq_index + 1 not in ptm_dict[ptm_name]['positions'].keys():
                    ptm_dict[ptm_name]['positions'][seq_index + 1] = []
                ptm_dict[ptm_name]['positions'][seq_index + 1].append(position)


# output - sequences numbered 0-based (and so are positions)
def find_ptm_sites(fasta_seqs_degapped, seq_ids):
    ptm_dict = {"ptm_phosph0": {'positions': {}},
                "ptm_phosph1": {'positions': {}},
                "ptm_phosph2": {'positions': {}},
                "ptm_phosph3": {'positions': {}},
                "ptm_phosphP": {'positions': {}},
                "ptm_acet0": {'positions': {}},
                "ptm_acet1": {'positions': {}},
                "ptm_acet2": {'positions': {}},
                "ptm_acet3": {'positions': {}},
                "ptm_Nglyc0": {'positions': {}},
                "ptm_Nglyc1": {'positions': {}},
                "ptm_Nglyc2": {'positions': {}},
                "ptm_Nglyc3": {'positions': {}},
                "ptm_amid0": {'positions': {}},
                "ptm_amid1": {'positions': {}},
                "ptm_amid2": {'positions': {}},
                "ptm_amid3": {'positions': {}},
                "ptm_hydroxy0": {'positions': {}},
                "ptm_hydroxy1": {'positions': {}},
                "ptm_hydroxy2": {'positions': {}},
                "ptm_hydroxy3": {'positions': {}},
                "ptm_methyl0": {'positions': {}},
                "ptm_methyl1": {'positions': {}},
                "ptm_methyl2": {'positions': {}},
                "ptm_methyl3": {'positions': {}},
                "ptm_Oglyc0": {'positions': {}},
                "ptm_Oglyc1": {'positions': {}},
                "ptm_Oglyc2": {'positions': {}},
                "ptm_Oglyc3": {'positions': {}},
                }
    for i, uni_id in enumerate(seq_ids):
        if uni_id:
            seq_data = get_uniprot_data(uni_id)
            parse_seq_data(seq_data, ptm_dict, i)
        predicted_phosph = run_netphos(fasta_seqs_degapped[(i * 2) + 1])
        ptm_dict["ptm_phosphP"]['positions'][i + 1] = predicted_phosph
    return ptm_dict


# output - sequences numbered 0-based (and so are positions)
def get_annotated_motifs(seq_ids):
    motif_dict = {}
    # get annotated motifs first
    for i, id_i in enumerate(seq_ids):
        if id_i:
            url = "http://elm.eu.org/instances.gff?q={}".format(id_i)
            try:
                req = urllib2.Request(url)
                response = urllib2.urlopen(req)
                features = response.read().splitlines()

                for j in features:
                    if 'sequence_feature' in j:
                        start = int(j.split()[3])
                        end = int(j.split()[4])
                        elm_id = 'motif_' + j.split()[8].split('=')[1]
                        if elm_id not in motif_dict:
                            motif_dict[elm_id] = {'prob': 1.0,
                                                  'positions': {}}
                        if i + 1 not in motif_dict[elm_id]['positions']:
                            motif_dict[elm_id]['positions'][i + 1] = []
                        motif_dict[elm_id]['positions'][i + 1].extend(
                            range(start, end + 1))
            except (urllib2.HTTPError,  urllib2.URLError,
                    SocketError, BadStatusLine):
                print 'get_annotated_motifs: HTTPError'
    return motif_dict


# find putative motifs by regular expressions
# TODO: implement filtering in structured regions
def get_predicted_motifs(sequences, elm_db, filter_motifs,):
    motif_dict = {}
    for i, seq_i in enumerate(sequences[1::2]):
        for elm_id in elm_db.keys():
            slim_j = elm_db[elm_id]
            reg = slim_j["comp_reg"]
            for match in reg.finditer(seq_i):
                m_sp = (match.span())
                prob = 1 + 1/math.log(slim_j["prob"], 10)
                start = m_sp[0] + 1
                end = m_sp[1] + 1
                if elm_id not in motif_dict:
                    motif_dict[elm_id] = {"prob": prob, "positions": {}}
                if i + 1 not in motif_dict[elm_id]["positions"]:
                    motif_dict[elm_id]["positions"][i + 1] = []
                motif_dict[elm_id]["positions"][i + 1].extend(range(start, end))
    return motif_dict


def find_motifs(sequences, seq_ids, elm_db):
    motif_dict = {}
    annotated = get_annotated_motifs(seq_ids)
    filter_motifs = False
    predicted = get_predicted_motifs(sequences, elm_db, filter_motifs)
    motif_dict = merge_motif_dicts(annotated, predicted)
    return motif_dict


def merge_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z


def find_domains(sequences, seq_ids):
    pass


def annotate(inname, outname):
    elm_db = parse_elm_db()
    fasta_seqs = read_fasta(inname)
    # check if it's aligned fasta
    if any(['-' in i for i in fasta_seqs[1::2]]):
        fasta_seqs_degapped = remove_gaps(fasta_seqs)
    else:
        fasta_seqs_degapped = fasta_seqs[:]
    seq_ids = find_uniprot_ids(fasta_seqs_degapped)
    ptm_sites = find_ptm_sites(fasta_seqs_degapped, seq_ids)
    motifs = find_motifs(fasta_seqs_degapped, seq_ids, elm_db)
    # domains = find_domains(fasta_seqs_degapped, seq_ids)
    metadata = merge_dicts(motifs, ptm_sites)
    print metadata
    write_conf_file(metadata, outname)


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('convert')
    parser = argparse.ArgumentParser(description='Convert fasta to a KMAD'
                                                 + ' compatible format')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')
    args = parser.parse_args()
    annotate(args.input_filename, args.output_filename)
