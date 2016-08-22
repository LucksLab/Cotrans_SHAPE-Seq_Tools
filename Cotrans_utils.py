"""
Modules used in creating reactivity and read matrices for Cotranscriptional SHAPE-Seq data.

Written by Angela M Yu <amy35@cornell.edu>, 2014-2016
Last edited: 5/4/2016
"""

import os
import glob
import re
import shutil


def format_rna_string(nt):
    return nt.upper().replace('T', 'U')


def end_match_strip(seq, adapter):
    """
    Iteratively checks if end of seq matches the adapter and removes the match
    from the end of seq. If match is not found, the last character in adapter
    is removed until a match is found or all characters in adapter are removed.
    """
    for i in range(0, len(adapter)):
        substr = adapter[:-i]
        if i == 0:
            substr = adapter
        if seq.endswith(substr):
            return seq[:-len(substr)], len(adapter) - i
    return seq, -len(seq)


def make_seq(seq, seqfilename):
    seq = format_rna_string(seq)

    with open(seqfilename, "w") as new_seq:
        line_to_write = ';\n{0}\n%s1'.format(seqfilename) % (seq)
        new_seq.write(line_to_write)

    return seqfilename


def remove_file(f):
    """ Removes a single file f """
    try:
        os.remove(f)
    except OSError:
        pass


def remove_files_with_ext(directory, extension):
    """ Remove all files in the specified directory with the extension """
    rmfiles = glob.glob(directory + "/*" + extension)
    for f in rmfiles:
        os.remove(f)


def recalc_thetas(thetas, start_i, end_i):
    """
    Recalculates theta values for a subsequence of the thetas.
    It assumes that the given sequence of thetas is already normalized.
    """
    cut_thetas = [float(r) for r in thetas[start_i:end_i]]
    cut_thetas_sum = sum([float(r) for r in cut_thetas])
    if cut_thetas_sum == 0:
        return cut_thetas
    recalced_thetas = [r/cut_thetas_sum for r in cut_thetas]
    return recalced_thetas


def calc_rho_from_theta_list(theta):
    rho = [float(t) * len(theta) for t in theta]
    return rho


def parse_reactivity_rho(infile, adapterseq, outputfile, endcut=0):
    """
    Parses reactivity file and outputs .theta and .seq, stripping off the
    adapter sequence (iteratively shrinking). Returns (positions, thetas,
    nt_sequence).
    """
    try:
        with open(infile, 'r') as f:
            f.readline()  # throw out header
            treated_sum, untreated_sum = [[int(s)] for s in re.split("\s+", f.readline())[4:6]]  # throw out nt = 0 except for rc
            seq = []
            theta = []
            theta_cut = []
            pos = []
            adapterseq = format_rna_string(adapterseq)

            for l in f:  # parse through file
                vars = re.split('\s+', l)
                seq.append(format_rna_string(vars[3]))
                theta.append(float(vars[7]))
                pos.append(vars[2])
                untreated_sum.append(int(vars[5]))
                treated_sum.append(int(vars[4]))

            seqstring = "".join(seq)
            (seq_cut, adapter_len) = end_match_strip(seqstring, adapterseq)

            # Also, endcut can be used to represent the polymerase footprint and remove these
            # bases from being used in the subsequent calculations.
            if endcut < 0:
                seq_cut = seq_cut[:endcut]
                adapter_len -= endcut
            pos = pos[:-adapter_len]
            untreated_sum = untreated_sum[:-adapter_len]
            treated_sum = treated_sum[:-adapter_len]

            theta_cut = [str(t) for t in recalc_thetas(theta, 0, -adapter_len)]
            rho_cut = calc_rho_from_theta_list(theta_cut)
            try:
                with open(outputfile+".theta", 'w') as out:
                    out.write("\n".join(["\t".join(z) for z in zip(pos, theta_cut)]))
            except EnvironmentError:
                print "Error opening output .theta file: " + outputfile + ".theta"

            try:
                with open(outputfile+".rho", 'w') as out:
                    out.write("\n".join(["\t".join([str(zi), str(zr)]) for zi,zr in zip(pos, rho_cut)]))
            except EnvironmentError:
                print "Error opening output .rho file: " + outputfile + ".rho"

            try:
                make_seq("".join(seq_cut), outputfile+".seq")
            except EnvironmentError:
                print "Error opening output .seq file: " + outputfile + ".seq"

            return (pos, rho_cut, theta)

    except EnvironmentError:
        print "Error opening reactivities file: " + infile


def reactivities_to_rho_file(input_dir, adapterseq, output_dir, rm_temp=True, min_len=0, max_len=0, endcut=0):
    """
    Takes a directory of reactivities files and the adapter sequence and
    creates .theta, .rho, and .seq files. Also outputs rho_table.txt
    which is a tab delimited file of rhos. Be aware that there is numerical
    error from using floats that are present in the output. There is an option (rm_temp)
    if you do not want to keep the .theta, .rho, and .seq files.
    """
    infiles = glob.glob(input_dir + "/*_reactivities.txt")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    outname = output_dir + "rho"

    rhos = {}
    for f in infiles:
        fname = re.findall("([^/]+).txt$", f)
        output_file_prefix = output_dir+"/"+fname[0]
        pos, rho_cut, theta = parse_reactivity_rho(f, adapterseq, output_file_prefix, endcut)
        if sum([float(t) for t in theta]) == 0:
            print("Not enough alignments in length {0} to calculate theta/rho. Line left blank.".format(len(pos)))
            rhos[len(pos)] = "-1\t"*len(pos) + "\n"
        else:
            rhos[len(pos)] = "\t".join([str(t) for t in rho_cut]) + "\n"

    outname += "_{0}min".format(min_len) if min_len != 0 else ""
    outname += "_{0}max".format(max_len) if max_len != 0 else ""
    outname += "_table.txt"
    with open(outname, 'w') as f:
        for key in sorted(rhos):
            if min_len == 0 or key >= min_len:    #Only apply logic if a value was supplied, but use bounds if so
                if max_len == 0 or key <= max_len:
                    f.write(rhos[key])
    if rm_temp:
        remove_file(input_dir + outname.split("/")[-1])
        shutil.move(outname, input_dir)
        remove_files_with_ext(output_dir, ".theta")
        remove_files_with_ext(output_dir, ".rho")
        remove_files_with_ext(output_dir, ".seq")
        remove_files_with_ext(output_dir, ".txt")
        os.rmdir(output_dir)
        return input_dir + os.path.basename(outname)
    return outname


def reactivities_to_reads_files(input_dir, adapterseq, min_len=0, max_len=0, endcut=0):
    """
    Takes a directory of reactivities files and the adapter sequence and
    creates tab delimited reads files:
    treated_mods_reads_table.txt, untreated_mods_reads_table.txt
    """
    infiles = glob.glob(input_dir + "/*_reactivities.txt")
    outname_pos = input_dir + "treated_mods_reads"
    outname_neg = input_dir + "untreated_mods_reads"
    outname_pos += "_{0}min".format(min_len) if min_len != 0 else ""
    outname_neg += "_{0}min".format(min_len) if min_len != 0 else ""
    outname_pos += "_{0}max".format(max_len) if max_len != 0 else ""
    outname_neg += "_{0}max".format(max_len) if max_len != 0 else ""
    outname_pos += "_table.txt"; outname_neg += "_table.txt";

    reads = {}
    for file in infiles:
        with open(file, 'r') as f:
            f.readline()
            lines = f.readlines()
            seq = "".join([a.split()[3] for a in lines[1:]])
            (seq_cut, adapter_len) = end_match_strip(seq, adapterseq)
            if endcut < 0:
                seq_cut = seq_cut[:endcut]
                adapter_len -= endcut

            if min_len == 0 or len(seq_cut) >= min_len:    #Only apply logic if a value was supplied, but use bounds if so
                if max_len == 0 or len(seq_cut) <= max_len:
                     reads[len(seq_cut)] = [a.split()[4:6] for a in lines[:-adapter_len]]
    with open(outname_pos, 'w') as f:
        with open(outname_neg, 'w') as f2:
            for key in sorted(reads):
                f.write("\t".join([a[0] for a in reads[key]]) + "\n")
                f2.write("\t".join([a[1] for a in reads[key]]) + "\n")

    return (outname_pos, outname_neg)
