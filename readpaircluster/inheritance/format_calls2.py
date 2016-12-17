#!/usr/bin/env python

"""
format_calls.py

Reformat SV caller output to bed or bedpe for cross-caller comparison.

Deletion, duplication, and inversion calls are reformatted as BED intervals,
and a separate BED file is output for each SV type. Translocation calls are
reformatted and output as a BEDPE file.

One-ended inversion and translocation breakpoints are output in separate
files for independent comparison.

Delly, lumpy-sv, gasvPro, softsearch, cnvnator, and the Talkowski lab
classifier are all supported call formats.
"""


import os
import argparse
import vcf
import gzip


# def vcf_to_bed(fin, prefix, outdir):
    # """
    # Reformat vcf output to bed or bedpe file

    # Deletions, duplications, and inversions are output to a bed
    # """
    # reader = vcf.Reader(fin)
    # try:
        # record = reader.next()
    # except ValueError as e:
        # # Improperly formatted VCF, replace 'NA' with '.'
        # if e.message == 'could not convert string to float: NA':
            # import sys
            # import traceback
            # traceback.print_exc()
            # sys.tracebacklimit = 0
            # raise ValueError('Improperly formatted VCF. This occurs when an '
                             # 'INFO field has been declared as a Float but one '
                             # 'or more records were assigned \'NA\' in this '
                             # 'field. Replace all such instances with \'.\' to '
                             # 'fix. `sed -i \'s/=NA/=./\'` is a quick fix.')

    # fin.seek(0)
    # reader = vcf.Reader(fin)

    # fouts = {}
    # for svtype in reader.alts.keys():
        # fouts[svtype] = open('%s/%s_%s.%s.bed' % (outdir, prefix, caller, svtype.lower()), 'w')

    # for record in reader:
        # svtype = record.INFO['SVTYPE']
        # if svtype in fouts.keys():
            # fout = fouts[svtype]
        # else:
            # continue

        # fout.write('%s\t%d\t%d\n' %
                   # (record.CHROM, record.POS, record.INFO['END']))

def vcf_to_bed(fin, prefix, outdir):
    """
    Reformat vcf output to bed or bedpe file

    Deletions, duplications, and inversions are output to a bed
    """
    reader = vcf.Reader(fin)
    try:
        record = reader.next()
    except ValueError as e:
        # Improperly formatted VCF, replace 'NA' with '.'
        if e.message == 'could not convert string to float: NA':
            import sys
            import traceback
            traceback.print_exc()
            sys.tracebacklimit = 0
            raise ValueError('Improperly formatted VCF. This occurs when an '
                             'INFO field has been declared as a Float but one '
                             'or more records were assigned \'NA\' in this '
                             'field. Replace all such instances with \'.\' to '
                             'fix. `sed -i \'s/=NA/=./\'` is a quick fix.')

    fin.seek(0)
    reader = vcf.Reader(fin)

    fouts = {}
    for svtype in reader.alts.keys():
        fouts[svtype] = open('%s/%s_%s.%s.bed' % (outdir, prefix, caller, svtype.lower()), 'w')

    for record in reader:
        svtype = record.INFO['SVTYPE']
        if svtype in fouts.keys():
            fout = fouts[svtype]
        else:
            continue

        fout.write('%s\t%d\t%d\n' %
                   (record.CHROM, record.POS, record.INFO['END']))


def delly_to_bed(**kwargs):
    """
    Reformat Delly output to bed or bedpe file

    Deletions, duplications, and inversions are output to a bed
    """
    reader = vcf.Reader(kwargs['fin'])
    prefix = kwargs['prefix']
    svtype = kwargs['svtype']
    outdir = kwargs['outdir']

    # if svtype != 'tloc' and svtype != 'inv':
    if False:
        # fout = open('%s/%s_delly.%s.bed' % (outdir, prefix, svtype), 'w')
        fout = open('%s/%s.%s.bed' % (outdir, prefix, svtype), 'w')
        for record in reader:
            fout.write("%s\t%d\t%d\t%s\n" % (record.CHROM, record.POS,
                                             record.INFO['END'], record.ID))
    else:
        fout = open('%s/%s.%s.bedpe' % (outdir, prefix, svtype), 'w')
        for record in reader:
            startA = record.POS + record.INFO['CIPOS'][0]
            endA = record.POS + record.INFO['CIPOS'][1]
            startB = record.INFO['END'] + record.INFO['CIEND'][0]
            endB = record.INFO['END'] + record.INFO['CIEND'][1]

            strandA, strandB = ['+' if strand == '5' else '-'
                                for strand in record.INFO['CT'].split('to')]

            fout.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\t.\t%s\t%s\n" %
                       (record.CHROM, startA, endA,
                        record.INFO['CHR2'], startB, endB,
                        record.ID, strandA, strandB))


def pindel_to_bed(**kwargs):
    """
    Reformat Delly output to bed or bedpe file

    Deletions, duplications, and inversions are output to a bed
    """
    reader = vcf.Reader(kwargs['fin'])
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    fouts = {}

    fouts['DEL'] = open('%s/%s_pindel.del.bed' % (outdir, prefix), 'w')
    fouts['DUP:TANDEM'] = open('%s/%s_pindel.dup.bed' % (outdir, prefix), 'w')
    fouts['INV'] = open('%s/%s_pindel.inv.bed' % (outdir, prefix), 'w')

    for record in reader:
        svtype = record.INFO['SVTYPE']
        if svtype in fouts.keys():
            fout = fouts[svtype]
        else:
            continue

        fout.write('%s\t%d\t%d\n' %
                   (record.CHROM, record.POS, record.INFO['END']))


def breakseq_to_bed(**kwargs):
    """
    Reformat Delly output to bed or bedpe file

    Deletions, duplications, and inversions are output to a bed
    """
    # reader = gzip.open(kwargs['fin'].name)
    reader = vcf.Reader(kwargs['fin'])
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    fouts = {}

    fouts['DEL'] = open('%s/%s_breakseq.del.bed' % (outdir, prefix), 'w')
    fouts['INS'] = open('%s/%s_breakseq.ins.bed' % (outdir, prefix), 'w')

    for record in reader:
        svtype = record.INFO['SVTYPE']
        if svtype in fouts.keys():
            fout = fouts[svtype]
        else:
            continue

        fout.write('%s\t%d\t%d\n' %
                   (record.CHROM, record.POS, record.INFO['END']))


def lumpy_new_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']
    svtype = kwargs['svtype']

    if svtype == 'del':
        fout = open('%s/%s_lumpy.del.bed' % (outdir, prefix), 'w')
    elif svtype == 'dup':
        fout = open('%s/%s_lumpy.dup.bed' % (outdir, prefix), 'w')
    else:
        raise Exception('Only dels and dups eligible for bed output')

    for line in fin:
        line = line.strip().split()
        chrom = line[0]
        min_start = int(line[1])
        min_end = int(line[4])

        fields = line[12].split(';')
        fields = dict([tuple(x.split('=')) for x in fields if len(x.split('=')) == 2])

        start_offset = abs(int(fields['CIPOS'].split(',')[0]))
        end_offset = abs(int(fields['CIEND'].split(',')[0]))

        start = min_start + start_offset + 1
        end = min_end + end_offset + 1

        ID = line[6]

        fout.write('%s\t%d\t%d\t%s\n' % (chrom, start, end, ID))


def lumpy_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    del_fout = open('%s/%s_lumpy.del.bed' % (outdir, prefix), 'w')
    dup_fout = open('%s/%s_lumpy.dup.bed' % (outdir, prefix), 'w')
    inv_fout = open('%s/%s_lumpy.inv.bed' % (outdir, prefix), 'w')
    inv_single = open('%s/%s_lumpy.inv.single.bed' % (outdir, prefix), 'w')
    tloc_fout = open('%s/%s_lumpy.tloc.bedpe' % (outdir, prefix), 'w')
    tloc_single = open('%s/%s_lumpy.tloc.single.bedpe' % (outdir, prefix), 'w')

    for line in fin:
        line = line.strip().split()

        chrom = line[0]
        ID = line[6]
        max_pos = line[13]
        start, end = [pos.split(':')[-1] for pos in max_pos.split(';')]

        if line[10] == 'TYPE:DELETION':
            del_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))
        elif line[10] == 'TYPE:DUPLICATION':
            dup_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))

        elif line[10] == 'TYPE:INVERSION':
            # If multiple strands exist, output to regular file
            if ';' in line[12]:
                inv_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))

            # Else output to the single ender output
            else:
                inv_single.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))

        elif line[10] == 'TYPE:INTERCHROM':
            # If multiple strands exist, output to regular file
            if ';' in line[12]:
                tloc_fout.write("\t".join(line) + "\n")

            # Else output to the single ender output
            else:
                tloc_single.write("\t".join(line) + "\n")
        else:
            print "Event %s not valid type (%s)" % (line[6], line[10])


def svsim_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    del_fout = open('%s/%s_svsim.del.bed' % (outdir, prefix), 'w')
    dup_fout = open('%s/%s_svsim.dup.bed' % (outdir, prefix), 'w')
    inv_fout = open('%s/%s_svsim.inv.bed' % (outdir, prefix), 'w')
    tloc_fout = open('%s/%s_svsim.tloc.bedpe' % (outdir, prefix), 'w')

    for line in fin:
        line = line.strip().split()

        chrom = line[0]
        ID = line[6]
        start = line[1]
        end = line[4]

        if 'DEL' in ID:
            del_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))
        elif 'DUP' in ID:
            dup_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))
        elif 'INV' in ID:
            inv_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end, ID))
        elif 'INR' in ID:
            tloc_fout.write("\t".join(line) + "\n")
        else:
            print "Event %s not valid type (%s)" % (line[6], line[10])


def cnvnator_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    del_fout = open('%s/%s_cnvnator.del.bed' % (outdir, prefix), 'w')
    dup_fout = open('%s/%s_cnvnator.dup.bed' % (outdir, prefix), 'w')

    line_num = 1
    for line in fin:
        line = line.strip().split()
        chrom = line[1].split(':')[0]
        start, end = line[1].split(':')[1].split('-')

        if line[0] == 'deletion':
            del_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end,
                                                 "\t".join(line[3:])))
        elif line[0] == 'duplication':
            dup_fout.write("%s\t%s\t%s\t%s\n" % (chrom, start, end,
                                                 "\t".join(line[3:])))
        else:
            print "Event on line %d not valid type. (%s)" % (line_num, line[0])

        line_num += 1


def classifier_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    svtype = kwargs['svtype']
    outdir = kwargs['outdir']

    if svtype != 'tloc':
        fout = open('%s/%s_classifier.%s.bed' % (outdir, prefix, svtype), 'w')

        # Inversions can have single enders
        if svtype is 'inv':
            fout_single = open('%s/%s_classifier.inv.single.bed' %
                               (outdir, prefix), 'w')

        for line in fin:
            line = line.strip().split()
            chrom = line[0]
            start = sum(int(x) for x in line[1:3]) / 2
            end = sum(int(x) for x in line[4:6]) / 2
            ID = line[6]

            if svtype is 'del' or 'dup':
                fout.write('%s\t%d\t%d\t%s\n' % (chrom, start, end, ID))

            # If it's an inversion, output single enders separately
            else:
                if line[8] is '.':
                    fout.write('%s\t%d\t%d\t%s\n' % (chrom, start, end, ID))
                else:
                    fout_single.write('%s\t%d\t%d\t%s\n' %
                                      (chrom, start, end, ID))

    else:
        fout = open('%s/%s_classifier.tloc.bedpe' % (outdir, prefix), 'w')
        fout_single = open('%s/%s_classifier.tloc.single.bedpe' %
                           (outdir, prefix), 'w')

        for line in fin:
            if line.split()[8] is '.':
                fout.write(line)
            else:
                fout_single.write(line)


def gasvpro_to_bed(**kwargs):
    fin = kwargs['fin']
    prefix = kwargs['prefix']
    outdir = kwargs['outdir']

    tlocs = ['TR+', 'TR-', 'TN+1', 'TN-1', 'TN+2', 'TN-2']
    invs = ['IR', 'I+', 'I-']
    tloc_fout = open('%s/%s_gasvpro.tloc.bedpe' % (outdir, prefix), 'w')
    inv_fout = open('%s/%s_gasvpro.inv.bedpe' % (outdir, prefix), 'w')
    tloc_single = open('%s/%s_gasvpro.tloc.single.bedpe' %
                       (outdir, prefix), 'w')

    fout = {
        'D': open('%s/%s_gasvpro.del.bed' % (outdir, prefix), 'w'),
        'IR': inv_fout,
        'I+': inv_fout,
        'I-': inv_fout,
        'TR+': tloc_fout,
        'TR-': tloc_fout,
        'TN+1': tloc_single,
        'TN-1': tloc_single,
        'TN+2': tloc_single,
        'TN-2': tloc_single
    }

    for line in fin:
        line = line.strip().split()

        # Skip "divergent" reads
        if line[7] == 'V':
            continue

        # If it's not a translocation, output as bed
        if line[7] not in tlocs and line[7] not in invs:
            chrom = line[1]
            start = sum([int(x) for x in line[2].split(',')]) / 2
            end = sum([int(x) for x in line[4].split(',')]) / 2
            ID = line[0]
            
            if chrom == '23':
                chrom = 'X'
            if chrom == '24':
                chrom = 'Y'

            fout[line[7]].write("%s\t%d\t%d\t%s\n" % (chrom, start, end, ID))

        # Else output as bedpe
        else:
            ID = line[0]
            chromA = line[1]
            startA, endA = line[2].split(',')
            chromB = line[3]
            startB, endB = line[4].split(',')
            size = line[5]
            
            if chromA == '23':
                chromA = 'X'
            if chromA == '24':
                chromA = 'Y'
            if chromB == '23':
                chromB = 'X'
            if chromB == '24':
                chromB = 'Y'

            svtype = line[7]
            if svtype.startswith('TR') or svtype == 'IR':
                strandA = '.'
                strandB = '.'
            elif svtype == 'TN+1':
                strandA = '+'
                strandB = '-'
            elif svtype == 'TN+2':
                strandA = '-'
                strandB = '+'
            elif svtype == 'TN-1' or svtype == 'I+':
                strandA = '+'
                strandB = '+'
            elif svtype == 'TN-2' or svtype == 'I-':
                strandA = '-'
                strandB = '-'

            fout[svtype].write(
                '%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t%s\n' %
                (chromA, startA, endA, chromB, startB, endB, ID,
                 strandA, strandB, size))


def softsearch_to_bed(args):
    pass


def seed_fout():
    pass


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'prog', choices=['delly', 'lumpy', 'gasvpro', 'cnvnator', 'softsearch',
                         'classifier', 'svsim', 'pindel', 'breakseq',
                         'genomestrip', 'lumpy_new'],
        metavar='prog',
        help="Program output to reformat. {delly,lumpy,gasvpro,cnvnator,"
        "classifier,svsim,pindel,breakseq,genomestrip,lumpy_new}")
    parser.add_argument('fin', type=argparse.FileType('r'),
                        help='SV caller output to reformat')
    parser.add_argument('prefix',
                        help="Prefix used in reformatted filenames (sample or "
                        "cohort name recommended). Output filenames will be "
                        "of the form prefix_prog.svtype.bed")

    parser.add_argument('-s', '--svtype', default=None,
                        choices=['del', 'dup', 'inv', 'tloc'],
                        metavar='TYPE',
                        help="Type of SV output to reformat. "
                        "Required for delly and classifier output. "
                        "{del,dup,inv,tloc}")
    parser.add_argument('-o', '--output-dir', default=os.getcwd(),
                        help='Output directory [CWD]')
    parser.add_argument('-f', '--format', default='bed',
                        choices=['bed','bedpe','vcf'],
                        help="File format to output. Default is bed for "
                        "del/dup and bedpe for inv/tloc. VCF not yet "
                        "supported.")

    args = parser.parse_args()

    if ((args.prog == 'delly' or args.prog == 'classifier' or args.prog == 'lumpy_new') and
            args.svtype is None):
        parser.error("Delly and classifier output require --svtype")

    formatters = {
        'delly': delly_to_bed,
        'lumpy': lumpy_to_bed,
        'classifier': classifier_to_bed,
        'gasvpro': gasvpro_to_bed,
        # 'softsearch': softsearch_to_bed,
        'cnvnator': cnvnator_to_bed,
        'svsim': svsim_to_bed,
        'pindel': pindel_to_bed,
        'breakseq': breakseq_to_bed,
        'lumpy_new': lumpy_new_to_bed
    }

    formatters[args.prog](fin=args.fin, prefix=args.prefix, svtype=args.svtype,
                          outdir=args.output_dir)

if __name__ == '__main__':
    main()
