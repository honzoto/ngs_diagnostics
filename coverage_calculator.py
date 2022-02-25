#!/home/sbsuser/miniconda3/bin/python

import sys, os
from pathlib import Path

# for data processing
from Bio import SeqIO
import numpy as np
import pysam

# for visualizations
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

"""
NOTES
v1  - added coverage / depth information to plots
    - aesthetic updates (colour, plot layout, x/y-axis font sizing)
    - detailed log file now includes bin location and sizes of each bin
    - frame correction for 1-based files
    - added average depth line to the plot

---
TO-DO / KNOWN ISSUES
- yet to develop ability to read subset files (i.e. if we don't want the whole genome coverage)
- if we knew the length of the reference beforehand, we shouldn't have to read the whole reference genome into file
- pysam module can only be used for linux, windows compatibility is in the works
- sam files are 1-based, whereas bam files are 0-based, make sure to correct for that

- coordinates file should be 0-based
"""


str_program = __file__.split("/")[-1]
str_version = "1.4"
str_help = """
[DESCRIPTION]
 This program accepts a SAM/BAM file as input and 

[USAGE]
 coverage_calculator.py -r <reference genome> -f <sam/bam file>

[HELP MENU]
 -h/--help          : shows this menu
 -f/--file          : path of the SAM/BAM file
 -r/--reference     : path of the reference FASTA file
 -c/--coords        : choose which coordinates to display (default=all)
 -b/--bins          : number of bins to show on the histogram (default=20)

"""


# ========================================[ GETTING ARGUMENTS ]=======================================

print("Welcome to honzo's {0} v{1} for SAM/BAM file visualization.".format(str_program, str_version))

# default arguments
pth_reference = Path("cov2_hg38.fasta")
pth_samfile = Path("COVID-0001_S1_L001_sars_aligned_sorted.bam")

pth_reference = None
pth_coords = "coord.csv"

str_pdfout = "plots_"+str(pth_samfile).replace(".bam", ".pdf")
tf_isbam = True
tf_useglobalcov = False # display coverage and depth for whole genome (true) or subset (false)
dic_params = {
    "bins": 50,
    "smoothing": "mean" # [sum/mean]
}

#plt.rcParams["font.family"] = "sans-serif"
#plt.rcParams.update({"font.sans-serif":"Arial"})
plt.rcParams["font.family"] = "Montserrat"

# RETRIEVE USER ARGUMENTS
if len(sys.argv) > 1:
    for i in range(len(sys.argv)):
        if sys.argv[i] == "-f" or sys.argv[i] == "--file":
            pth_samfile = Path(sys.argv[i+1])
            str_ext = pth_samfile.suffix
            if str_ext == ".sam":
                tf_isbam = False
        elif sys.argv[i] == "-r" or sys.argv[i] == "--reference":
            pth_reference = Path(sys.argv[i+1])
            pth_coords = None
        elif sys.argv[i] == "-c" or sys.argv[i] == "--coords":
            pth_coords = Path(sys.argv[i+1])
            pth_reference = None

        # optional arguments
        elif sys.argv[i] == "-b" or sys.argv[i] == "--bins":
            dic_params["bins"] = int(sys.argv[i+1])

        elif sys.argv[i] == "-h" or sys.argv[i] == "--help":
            print(str_help)
            quit()

str_basename = pth_samfile.stem
wf_log = open("log_coverage-calc.txt", "wt")

# TITLE BLOCK
wf_log.write("This is honzo's {0} for coverage visualization.\n".format(str_program))
wf_log.write("\nParameters:\n")
wf_log.write("(-f) Input file     : {0}\n".format(pth_samfile))
wf_log.write("(-r) Reference file : {0}\n".format(pth_reference))
wf_log.write("(-s) Subset file    : {0}\n".format(pth_coords))
wf_log.write("\n")

# =====================================[ GETTING PARAMS DETAILS ]=====================================

# getting reference lengths from file
dic_coords = {} # [accession] = [start, end]
dic_covmap = {} # a list of zero-values spanning the length of each sequence

if pth_reference: # option 1: build coverage map from supplied fasta file
    print("Getting length of references from file: {0}".format(pth_reference))
    wf_log.write("\nBuilding coverage map from reference file.\n")
    for record in SeqIO.parse(pth_reference, "fasta"):
        dic_covmap[record.id] = np.zeros(len(record.seq))
        print("Record: {0}\tSequence length: {1}".format(record.id, len(record.seq)))
        wf_log.write("Record: {0}, length: {1}\n".format(record.id, len(record.seq)))

elif pth_coords: # option 2: build coverage map from set coordinates in csv file
    with open(pth_coords, "rt") as rf_coords:
        print("Building coverage map from preset coordinates.")
        wf_log.write("\nBuilding coverage map from preset coordinates.\n")
        rf_coords.readline() # header line
        for line in rf_coords:
            lst_vals = line.strip().split(",")
            # dic_coords[accession] = (start, end)
            dic_coords[lst_vals[0]] = int(lst_vals[1]), int(lst_vals[2])
            dic_covmap[lst_vals[0]] = np.zeros(int(lst_vals[3]))
            wf_log.write("Record: {0}, start: {1}, end: {2}, genome length: {3}\n" \
                .format(lst_vals[0], lst_vals[1], lst_vals[2], len(dic_covmap[lst_vals[0]])))

else:
    print(str_help)
    print("[error] no reference file (-r) or coordinates file (-c) indicated.")


# =================================[ CALCULATE COVERAGE FROM BAM/SAM ]================================

# Indexing bam/sam file if necessary
if not os.path.exists(str(pth_samfile)+".bai"):
    print("\nCould not locate .bai/.sai file. Indexing file now.")
    wf_log.write("[warning] could not find .bai/.sai file. Sorting and indexing file.\n")

    str_sortedsam = str_basename+"_sorted"+str_ext
    cmd = "samtools sort -o {o} {i}".format(i=pth_samfile, o=str_sortedsam)
    wf_log.write("[cmd] {0}\n".format(cmd))
    print("[cmd]", cmd)
    os.system(cmd)

    cmd = "samtools index {i}".format(i=str_sortedsam)
    wf_log.write("[cmd] {0}\n".format(cmd))
    print("[cmd]", cmd)
    os.system(cmd)

    pth_samfile = str_sortedsam

rf_sam = pysam.AlignmentFile(pth_samfile, "rb")
print("Parsing samfile.")

for read in rf_sam.fetch():
    if not read.reference_name in dic_covmap: continue
    # we're enumerating, but do we want to get the number of entries in the samfile to show progress?
    # https://pysam.readthedocs.io/en/latest/usage.html
    start = read.reference_start - 1
    end = start + len(read.query_sequence)
    # +2 because (1) coords file is 1-based, final x-axis is 1-based
    if pth_coords:
        dic_covmap[read.reference_name][start+2:end+2] += 1
    else:
        dic_covmap[read.reference_name][start+1:end+1] += 1

# ====================================[ GENERATING VISUALIZATIONS ]===================================

with PdfPages(str_pdfout) as pdf:
    # each accession will have its own coverage map
    for str_accession in dic_covmap:
        print("Generating visualizations for accession: {0}".format(str_accession))
        wf_log.write("Accession: {0}\n".format(str_accession))

        arr_covmap = dic_covmap[str_accession]
        int_start = 0
        int_end = len(arr_covmap)
        lst_xvals = [] # for graphing
        lst_yvals = []

        if pth_coords:
            if str_accession in dic_coords:
                int_start = dic_coords[str_accession][0]
                int_end = dic_coords[str_accession][1]
            else:
                print("[warning] {0} not found in coordinates file.".format(str_accession))
                continue

        # in case we have more bins than we have # bases, we want to set the bins = # bases
        int_blen = (int_end - int_start) // dic_params["bins"]
        if int_blen == 0: 
            int_blen = 1
        
        print("Generating bins for [{0}], start={1}, end={2}".format(str_accession, int_start, int_end))
        for i in range(int_start, int_end, int_blen):
            lst_xvals.append(str(i))
            if dic_params["smoothing"] == "sum":
                lst_yvals.append(np.sum(arr_covmap[i:i+int_blen]))
            elif dic_params["smoothing"] == "mean":
                flt_mean = np.mean(arr_covmap[i:i+int_blen])
                lst_yvals.append(flt_mean)

        if tf_useglobalcov:
            flt_coverage = round(np.count_nonzero(arr_covmap) / len(arr_covmap) * 100, 2)
            flt_depth = round(np.mean(arr_covmap), 2)
        else:
            arr_yvals = np.array(lst_yvals)
            flt_coverage = round(np.count_nonzero(arr_yvals) / len(arr_yvals) * 100, 2)
            flt_depth = round(np.mean(arr_yvals), 2)

        wf_log.write("xvals: "+str(lst_xvals)+"\n")
        wf_log.write("yvals: "+str(lst_yvals)+"\n")
        wf_log.write("Coverage: {0}\n".format(flt_coverage))
        wf_log.write("Average read depth: {0}\n\n".format(flt_depth))

        # plotting
        fig, ax = plt.subplots(figsize= (14.0, 8.5))
        ax.bar(lst_xvals, lst_yvals, zorder=2, align='edge', width=0.95, color="#407BFF", alpha=0.8)
        ax.axhline(y=flt_depth, zorder=3, color="#FF6066", alpha=0.6)
        ax.set_title("{a} Coverage Map\n(est. coverage = {c}%, avg. depth={d})" \
            .format(a=str_accession, c=flt_coverage, d=flt_depth), fontsize=18)
        ax.set_ylabel("{0} counts per bin".format(dic_params["smoothing"]), fontsize=14)
        ax.set_xlabel("Genome region (bp)", fontsize=14)
        ax.tick_params(axis='x', labelrotation=90)
        if pth_coords:
            # if we create a subset, the x-axis does not cotnain the last tick
            ax.set_xticks(lst_xvals+[str(int_end)])

        # aesthetics
        ax.grid(which='major', color='#E2E2E2', linestyle='--', zorder=0)
        ax.set_facecolor('#F5F5F5')
        lst_sides = ["top","right"]
        for side in lst_sides:
            ax.spines[side].set_visible(False)

        plt.tight_layout()
        pdf.savefig()
        plt.close()

wf_log.close()
print("program complete.")