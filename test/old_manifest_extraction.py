# -*- coding: UTF-8 -*-
from __future__ import print_function
import sys
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, NormalizationTransform, code2genotype
from collections import namedtuple


def get_manifest(manifest_path, extraction_path, sep=','):
    """
    Extract manifest data about SNPs and write it to a ../data/extracted folder
    :param manifest_path: str - path to file
    :param extraction_path: str - path to directory where extracted files will be stored
    :param sep: str - separator which is used for writing, comma by default
    :return:
    """
    # Add {} to use it later in formatting names
    path_to_save = extraction_path + '/{}'
    # Get data
    manifest = BeadPoolManifest(manifest_path)

    # Check whether name of manifest is the same as coded in the file
    assert manifest.manifest_name == manifest_path.split('/')[-1], \
        "Name of manifest file doesn't match with manifest name from file"

    # List of fields which should be extracted from the manifest
    manifest_extract = ['names', 'chroms', 'map_infos', 'ref_strands', 'source_strands', 'snps']

    content = []
    # Iterate over attributes of manifest object
    for attr in manifest_extract:
        content.append((attr, map(str, manifest.__getattribute__(attr))))

    # Initialize variables
    name = manifest.manifest_name.split('.')[0] + '_old.csv'
    length = len(content[0][1])
    rows = []

    # Make header
    header = sep.join([content[i][0] for i in range(len(content))])

    # Create normal df structure
    for i in range(length):
        row = sep.join([content[j][1][i] for j in range(len(content))])
        rows.append(row)

    # Write to a file
    with open(path_to_save.format(name), 'w') as dest:
        dest.write(header + '\n' + '\n'.join(rows))


if __name__ == '__main__':
    # manifest_path = "/home/arleg/repos/cowsdb/data/manifest/QCArray_96XT_20005626_A1.bpm"
    # manifest_path = "/home/arleg/repos/bead_array_files/data/OvineSNP50v2_XT_20006795X356271_A1.bpm"
    manifest_path = sys.argv[1]
    get_manifest(manifest_path, extraction_path='/home/arleg/repos/bead_array_files/extracted')

