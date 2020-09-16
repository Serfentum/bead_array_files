import sys
# import argparse
import os
from datetime import datetime
from GenotypeCalls import *
from BeadPoolManifest import *


delim = "\t"
output_file = '/home/arleg/repos/bead_array_files/extracted/report1'
manifest_path = '/home/arleg/repos/bead_array_files/data/OvineSNP50v2_XT_20006795X356271_A1.bpm'
manifest= BeadPoolManifest(manifest_path)
gtc_directory = '/home/arleg/repos/bead_array_files/data'

with open(output_file, "w") as output_handle:
    output_handle.write("[Header]\n")
    output_handle.write(delim.join(["Processing Date", datetime.now().strftime("%m/%d/%Y %I:%M %p")])+ "\n")
    output_handle.write(delim.join(["Content", os.path.basename(manifest_path)]) + "\n")
    output_handle.write(delim.join(["Num SNPs", str(len(manifest.names))]) + "\n")
    output_handle.write(delim.join(["Total SNPs", str(len(manifest.names))]) + "\n")

    samples = []
    for gtc_file in os.listdir(gtc_directory):
        if gtc_file.lower().endswith(".gtc"):
            samples.append(gtc_file)

    output_handle.write(delim.join(["Num Samples", str(len(samples))]) + "\n")
    output_handle.write(delim.join(["Total Samples", str(len(samples))]) + "\n")

    output_handle.write("[Data]\n")
    output_handle.write(delim.join(["SNP Name", "Sample ID", "Chr", "MapInfo", "Alleles - AB", "Alleles - Plus", "Alleles - Forward"]) + "\n")
    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(gtc_directory, gtc_file)
        gtc = GenotypeCalls(gtc_file)
        genotypes = gtc.get_genotypes()
        plus_strand_genotypes = gtc.get_base_calls_plus_strand(manifest.snps, manifest.ref_strands)
        forward_strand_genotypes = gtc.get_base_calls_forward_strand(manifest.snps, manifest.source_strands)

        assert len(genotypes) == len(manifest.names)
        for (name, chrom, map_info, genotype, ref_strand_genotype, source_strand_genotype) in zip(manifest.names, manifest.chroms, manifest.map_infos, genotypes, plus_strand_genotypes, forward_strand_genotypes):
            output_handle.write(delim.join([name, os.path.basename(gtc_file)[:-4], chrom, str(map_info), code2genotype[genotype], ref_strand_genotype, source_strand_genotype]) + "\n")
