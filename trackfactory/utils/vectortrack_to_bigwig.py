'''
Created on Mar 19, 2011

@author: mkiyer
'''
import argparse
import logging

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", dest="ucsc_bin", default="bedGraphToBigWig")
    parser.add_argument("track", dest="file")
    parser.add_argument("file", dest="file")
    
    
    
    options = parser.parse_args()
    
    
if __name__ == '__main__':
    main()
