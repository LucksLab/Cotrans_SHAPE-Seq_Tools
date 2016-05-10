#!/usr/bin/env python
# -*- coding: utf-8 -*-

##Written by Kyle Watters on October, 7th 2014 to parse reactivities data into a matrix format for matlab
#Edited by Angela M Yu on May 4, 2016

import getopt
import sys
import Cotrans_utils

name = "Cotrans_plot_Matlab.py"

help_message='''

{0} takes a directory containing cotranscriptional SHAPE-Seq reactivities and converts it to a nice matrix form

Usage:                                                                                                                                                                                                                   
   python {0} [options] <reactivities_dir>                                                                                                                                              

General Options:                                                                                                                                                                                                        
-h, --help                  Opens help message
-v, --version               Displays version number
-l, --linker <string>       Linker sequence to exclude from RNAs (default is IDT mod: CTGACTCGGGCACCAAGGA)
-r, --recalc_end            Recalculate rhos. Specify the last index of the rhos to be considered. Ex. -1 will drop the last rho value and renormalize based off of reduced length.
-i, --min-len <N>           Minimum length to include (default is first length)
-m, --max-len <N>           Maximum length to include (default is last length)

'''.format(name)

def get_version():
    return "0.0.2"

class Usage(Exception):
    def __init__(self,msg):
        self.msg = msg

class Params:
    def __init__(self):
        pass

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hvl:r:i:m:",
                                       ["version",
                                        "help",
					"linker=",
                                        "min-len=",
                                        "max-len="])

        except getopt.error, msg:
            raise Usage(msg)
          
        linker_seq = "CTGACTCGGGCACCAAGGA"     
        recalc_end = 0
        min_len = 0
        max_len = 0
        for option, value in opts:
            if option in ("-v", "--version"):
                print "%s v%s" % (name,get_version())
                exit(0)
            elif option in ("-h", "--help"):
                raise Usage(help_message)
            elif option in ("-l","--linker"):
                linker_seq = value
            elif option in ("-i","--min-len"):
                min_len = int(value)
            elif option in ("-m","--max-len"):
                max_len = int(value)

	if len(args) != 1:
            raise Usage(help_message)
        
        return args,linker_seq,recalc_end,min_len,max_len


    def check(self):
        pass


def main(argv=None,):

    params = Params()

    try:
        if argv is None:
            argv = sys.argv
            
            args,linker_seq,recalc_end,min_len,max_len = params.parse_options(argv)
            params.check()
        reactivities_dir = args[0]
        output_dir = reactivities_dir + '/tmp_rho_conversion/'
        rm_temp = True
        Cotrans_utils.reactivities_to_rho_file(reactivities_dir, linker_seq, output_dir, rm_temp, min_len, max_len)
        Cotrans_utils.reactivities_to_reads_files(reactivities_dir, linker_seq, min_len, max_len)
        
    except Usage, err:	
	print sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg) 
	print >> sys.stderr, "" 
        return 2

if __name__ == "__main__":
    sys.exit(main())


