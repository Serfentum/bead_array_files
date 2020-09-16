import subprocess as sp
import shlex
import sys
import os
from new_manifest_extraction import get_manifest


def compare_bpm(manifest_path, extraction_path='/home/arleg/repos/bead_array_files/extracted'):
    """
    Check whether manifests obtained with 2 procedures are identical
    :param manifest_path: str - path to bpm
    :param extraction_path: str - path to dir, where extracted files should be stored
    :return: bool - True if files are identical
    """
    bpm = manifest_path.split('/')[-1]
    bpm = bpm.split('.')[0]

    sp.run(shlex.split(f'python2 /home/arleg/repos/bead_array_files/test/old_manifest_extraction.py {manifest_path}'))
    get_manifest(manifest_path, extraction_path)

    res = sp.run(shlex.split(
        f'diff -s {os.path.join(extraction_path, f"{bpm}_old.csv")} {os.path.join(extraction_path, f"{bpm}_new.csv")}'),
                 stdout=sp.PIPE)
    if res.stdout.decode().strip().endswith('are identical'):
        return True
    print(res.stdout.decode())
    return False


if __name__ == '__main__':
    bpms = ("/home/arleg/repos/bead_array_files/data/OvineSNP50v2_XT_20006795X356271_A1.bpm",
            "/home/arleg/repos/bead_array_files/data/HumanKaryomap-12v1_A.bpm")

    truthiness = []

    for bpm in bpms:
        is_passed = compare_bpm(manifest_path=bpm,
                                extraction_path='/home/arleg/repos/bead_array_files/extracted')
        print(bpm, is_passed, sep='\t')
        truthiness.append(is_passed)

        assert all(truthiness), 'Some of the cases failed'
