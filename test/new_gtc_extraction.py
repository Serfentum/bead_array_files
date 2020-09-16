from GenotypeCalls import *
from BeadPoolManifest import *
# NormalizationTransform, SourceStrand, code2genotype
from collections import namedtuple
import sys


# todo add assert with manifest name manifest.name and checksum of manifest and man in db - probably should be in main
#  py3 file
def extract(gtc_path, extraction_path, manifest_path="/home/ailin/repo_new/data/BovineSNP50_v3_A1.bpm"):
    """
    Extract genotyping data -
    ballele_freqs, base_calls, genotypes, genotype_scores, logr_ratios, raw_x_intensities, raw_y_intensities,
    normalized_intensities, names, chroms, map_infos, ref_strands, source_strands and snps -
    and write it to a file
    Also extract general sample information -
    call_rate, cluster_file, gender, imaging_date, autocall_date, scanner_data, snp_manifest, is_write_complete,
    sample_name, sample_plate, sample_well -
    and write it to .sinfo file
    :param gtc_path: str - path to gtc with genotyping data
    :param extraction_path: str - path to directory where extracted files will be stored
    :param manifest_path: str - path to manifest file used for creation of this gtc
    :return:
    """
    # Add {} to use it later in formatting names
    path_to_save = extraction_path + '/{}'
    # Get gtc and manifest objects
    gtc = GenotypeCalls(gtc_path)
    manifest = BeadPoolManifest(manifest_path)

    # Structure for ordered names and methods of gtc fields
    field = namedtuple('field', ['name', 'method'])

    # List of fields which should be extracted from gtc
    gtc_extract = [field('ballele_freqs', GenotypeCalls.get_ballele_freqs),
                   field('genotypes', GenotypeCalls.get_genotypes),
                   field('genotype_scores', GenotypeCalls.get_genotype_scores),
                   field('logr_ratios', GenotypeCalls.get_logr_ratios),
                   field('raw_x_intensities', GenotypeCalls.get_raw_x_intensities),
                   field('raw_y_intensities', GenotypeCalls.get_raw_y_intensities),
                   field('normalized_intensities', lambda x: GenotypeCalls.get_normalized_intensities(x, manifest.normalization_lookups))]

    # I don't see place in db where we use this data
    # List of fields which correspond to a whole sample and extracted from gtc
    sample_info = [field('call_rate', GenotypeCalls.get_call_rate),
                   field('cluster_file', GenotypeCalls.get_cluster_file),
                   field('gender', GenotypeCalls.get_gender),
                   field('imaging_date', GenotypeCalls.get_imaging_date),
                   field('autocall_date', GenotypeCalls.get_autocall_date),
                   field('scanner_data', GenotypeCalls.get_scanner_data),
                   field('snp_manifest', GenotypeCalls.get_snp_manifest),
                   field('is_write_complete', GenotypeCalls.is_write_complete),
                   field('sample_name', GenotypeCalls.get_sample_name),
                   field('sample_plate', GenotypeCalls.get_sample_plate),
                   field('sample_well', GenotypeCalls.get_sample_well)]

    # Containers for data
    content = []
    general_info = []

    # Get content from gtc
    # Iterate over fields which should be extracted in gtc, transform them to str
    for name, method in gtc_extract:
        res = method(gtc)
        # For normalized_intensities divide the list of (x, y) intensities into lists of x and y
        if name != 'normalized_intensities':
            if not isinstance(res, str):
                try:
                    res = list(map(str, res))
                except TypeError:
                    res = str(res)
            content.append((name, res))
        else:
            # Compute r and theta values
            polar = list(map(NormalizationTransform.rect_to_polar, res))
            content.append(('normalized_x_intensities', [str(x) for x, y in res]))
            content.append(('normalized_y_intensities', [str(y) for x, y in res]))
            content.append(('r', [str(r) for r, theta in polar]))
            content.append(('theta', [str(theta) for r, theta in polar]))

    # Get base calls and their forward encoding
    base_calls = GenotypeCalls.get_base_calls(gtc)
    # print(SourceStrand.Forward)
    genotype_forward = GenotypeCalls.get_base_calls_forward_strand(gtc,
                                                          base_calls,
                                                          [SourceStrand.Forward for i in range(len(base_calls))])
    # Write them to collection
    content.append(('base_calls', base_calls))
    content.append(('genotype_forward', genotype_forward))

    # Iterate over sample information attributes of gtc object
    for name, method in sample_info:
        res = str(method(gtc))
        general_info.append((name, res))


    # Initialize variables
    length = len(content[0][1])
    sep = ','
    rows = []



    # Make header
    header = sep.join([content[i][0] for i in range(len(content))])

    # try:
    # Create normal df structure
    for i in range(length):
        row = sep.join([content[j][1][i] for j in range(len(content))])
        rows.append(row)

    # File names
    name = gtc_path.split('/')[-1].split('.')[0]
    sinfo_name = name + '_new.sinfo'
    data_name = name + '_new.csv'

    # Write data to a file
    with open(path_to_save.format(data_name), 'w') as dest:
        dest.write(header + '\n' + '\n'.join(rows))

    # Write sample information to a file
    with open(path_to_save.format(sinfo_name), 'w') as dest:
        dest.write('\n'.join([sep.join(record) for record in general_info]))


if __name__ == "__main__":
    extract(gtc_path="/home/arleg/repos/bead_array_files/data/203470490019_R07C01 (2).gtc",
            extraction_path='/home/arleg/repos/bead_array_files/extracted',
            manifest_path='/home/arleg/repos/bead_array_files/data/OvineSNP50v2_XT_20006795X356271_A1.bpm')
