#!/usr/bin/env python3

import argparse
import gzip

NODE_FILE = './taxonomy/nodes_lite.dmp.gz'


def __parse_taxonomy_nodes():
    nodes = {}
    with gzip.open(NODE_FILE, 'rt') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 2:
                key = columns[0]
                value = columns[1]
                nodes[key] = value
    return nodes

if __name__ == '__main__':
    # I want to take as parameters a gzip blast results file, a donor taxid and a group taxid
    parser = argparse.ArgumentParser(description='Calculate the AI and AHS for a group of sequences')
    parser.add_argument('--blast_out', help='The blast results file')
    parser.add_argument('--donor_taxid', help='Ingroup taxid (default Metazoan)', default=33208, type=int)
    parser.add_argument('--group_taxid', help='Excluded group (default Rotifera)', default=10190, type=int)
    parser.add_argument('--output', help='Output file', default='alienness_results.tsv')
    args = parser.parse_args()
    #print(args.blast_results)
    print('--- Parsing taxonomy nodes...')
    nodes = __parse_taxonomy_nodes()
    print('Done')

    
    