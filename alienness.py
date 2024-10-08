#!/usr/bin/env python3

import argparse
import gzip
from math import exp, log

NODE_FILE = './taxonomy/nodes_lite.dmp.gz'
SIGNIFICANT_IDS_FILE = './taxonomy/significant_ids.gz'


def __parse_taxonomy_nodes():
    """
    Parse the NCBI node lite file to extract a dict of key: tax_id, value: parent_tax_id
    """
    nodes = {}
    with gzip.open(NODE_FILE, 'rt') as file:
        for line in file:
            lsplt = line.strip().split('\t')
            if len(lsplt) == 2:
                nodes[int(lsplt[0])] = int(lsplt[1])
    return nodes


def __parse_significant_ids():
    """
    Parse the significant ids
    """
    significant_ids = set()
    with gzip.open(SIGNIFICANT_IDS_FILE, 'rt') as file:
        for line in file:
            significant_ids.add(int(line.strip()))
    return significant_ids


def __parse_blast_results(blast_file):
    """
    Parse the BLAST result file to extract the list of e-value and bit-score with their associated taxid
    """
    blast_results = {}
    with gzip.open(blast_file, 'rt') as file:
        for line in file:
            lsplt = line.strip().split('\t')
            if len(lsplt) != 13 or lsplt[-1] == 'NA':
                continue
            if lsplt[0] not in blast_results:
                blast_results[lsplt[0]] = []
            blast_results[lsplt[0]].append((float(lsplt[10]), float(lsplt[11]), int(lsplt[12])))
    return blast_results


def is_in_group(nodes, taxid, group_taxid):
    """
    Return True if the taxid belongs to the group_taxid, False otherwise.
    """
    if taxid == group_taxid:
        return True
    elif taxid not in nodes or nodes[taxid] == taxid:
        return False
    else:
        return is_in_group(nodes, nodes[taxid], group_taxid)


def __compute_ai(blast_results, nodes, donor_taxid, excluded_taxid, significant_ids):
    """
    Compute the Alienness Index (AI) for a given sequence results according to https://alienness.sophia.inrae.fr/cgi/faq.cgi
    """
    best_donor, best_recipient = 1, 1
    for (evalue, _, taxid) in blast_results: 
        if taxid not in significant_ids or is_in_group(nodes, taxid, excluded_taxid):
            continue
        elif is_in_group(nodes, taxid, donor_taxid):
            best_recipient = min(best_recipient, evalue)      
        elif not is_in_group(nodes, taxid, donor_taxid):
            best_donor = min(best_donor, evalue)
    #print(f'Best donor: {best_donor}, Best recipient: {best_recipient}')
    if best_donor == 1 and best_recipient == 1:
        return 'NA'
    return round(log(best_recipient + exp(-200)) - log(best_donor + exp(-200)), 2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate the AI and AHS for a group of sequences')
    parser.add_argument('--blast_out', help='The blast results file')
    parser.add_argument('--donor_taxid', help='Ingroup taxid (default Metazoan)', default=33208, type=int)
    parser.add_argument('--excluded_taxid', help='Excluded group (default Rotifera)', default=10190, type=int)
    parser.add_argument('--output', help='Output file', default='alienness_results.tsv')
    args = parser.parse_args()
    print('--- Parsing taxonomy nodes...')
    nodes = __parse_taxonomy_nodes()
    print('--- Parsing significant ids...')
    significant_ids = __parse_significant_ids()
    print('--- Parsing BLAST results...')
    blast_results = __parse_blast_results(args.blast_out)
    print('--- Computing AI...')
    for protid in blast_results:
        print(f'{protid}: {__compute_ai(blast_results[protid], nodes, args.donor_taxid, args.excluded_taxid, significant_ids)}')
    print('Done')