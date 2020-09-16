import subprocess as sp
import shlex
import sys
import os
import numpy as np
import pandas as pd
from new_gtc_extraction import extract


def compare_gtc(gtc_path, manifest_path, extraction_path='/home/arleg/repos/bead_array_files/extracted'):
    """
    Check whether manifests obtained with 2 procedures are identical
    :param gtc_path: str - path to gtc
    :param manifest_path: str - path to bpm
    :param extraction_path: str - path to dir, where extracted files should be stored
    :return: bool - True if files are identical
    """
    gtc = gtc_path.split('/')[-1]
    gtc = gtc.split('.')[0]
    sp.run(shlex.split(f'python2 /home/arleg/repos/bead_array_files/test/old_gtc_extraction.py "{manifest_path}" "{gtc_path}"'))

    extract(gtc_path, extraction_path, manifest_path)

    # Read extracted data
    # os.path.join(extraction_path, f"{gtc}_old.csv")}
    old = pd.read_csv(os.path.join(extraction_path, f'{gtc}_old.csv'))
    new = pd.read_csv(os.path.join(extraction_path, f'{gtc}_new.csv'))

    old_numerical = old.select_dtypes(include=np.number)
    new_numerical = new.select_dtypes(include=np.number)

    old_qualitative = old.select_dtypes(exclude=np.number)
    new_qualitative = new.select_dtypes(exclude=np.number)

    # Whether numerical values are approximately equal
    num_equality = np.allclose(old_numerical.values, new_numerical.values, equal_nan=True)
    quality_equality = old_qualitative.equals(new_qualitative)

    if num_equality and quality_equality:
        return True
    return False


if __name__ == '__main__':
    bpm = "/home/arleg/repos/bead_array_files/data/OvineSNP50v2_XT_20006795X356271_A1.bpm"
    gtcs = ("/home/arleg/repos/bead_array_files/data/203470490019_R07C01 (2).gtc",
            "/home/arleg/repos/bead_array_files/data/203470490019_R07C02 (2).gtc",
            "/home/arleg/repos/bead_array_files/data/203470490019_R07C03 (2).gtc",
            "/home/arleg/repos/bead_array_files/data/203470490019_R07C04 (2).gtc")

    truthiness = []

    for gtc in gtcs:
        is_passed = compare_gtc(gtc_path=gtc,
                                manifest_path=bpm,
                                extraction_path='/home/arleg/repos/bead_array_files/extracted')
        print(gtc, is_passed, sep='\t')
        truthiness.append(is_passed)

        assert all(truthiness), 'Some of the cases failed'
