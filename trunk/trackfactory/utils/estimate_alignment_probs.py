'''
Created on Mar 7, 2011

@author: mkiyer
'''
import argparse
import logging
import os
import subprocess
import numpy as np
import pysam

DESCRIPTION=("estimate alignment probability of multimapping reads")

def estimate_alignment_probs(file):
    samfh = pysam.Samfile(file, "rb")
    # file must be name sorted
    samfh.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument("file")

if __name__ == '__main__':
    main()