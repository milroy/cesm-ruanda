#!/bin/python

import os
import sys
import glob
from bs4 import BeautifulSoup as bs
import copy
from collections import defaultdict
from argparse import ArgumentParser
import pdb

parser = ArgumentParser()
parser.add_argument("--cov", dest="covPath", required=True, type=str)
parser.add_argument("--src", dest="srcPath", required=True, type=str)

args = parser.parse_args()
covPath = args.covPath
srcPath = args.srcPath

def readF(path):
    with open(path, 'r') as f:
        for line in f:
            yield line

def commentSubprog(witer, line, output, spType, name, starts):
    innerStarts = copy.deepcopy(starts)
    endNameTarget = 'end' + ' ' + spType #+ ' ' + name
    # You can write endsubroutine!
    endNameTarget1 = 'end' + spType #+ ' ' + name
    output.write('! Subprogram not used ' + line)
    new = next(witer, None)
    endName = ''
    while new is not None and endNameTarget not in endName and endNameTarget1 not in endName:
        # This breaks for two line subprograms, so need to check this.
        if 'end module' not in ' '.join(new.strip().lower().split()) or 'endmodule' not in ' '.join(new.strip().lower().split()):
            output.write('! Subprogram not used ' + new)

        else:
            output.write(new)

        new = next(witer, None)

        if new:
            newProcessed = ' '.join(new.strip().lower().split())
            newSplits = newProcessed.split()

            if len(newSplits) > 2:
                endName = newSplits[0] + ' ' + newSplits[1] + ' ' + newSplits[2]
                if endNameTarget in endName or endNameTarget1 in endName or 'end module' in endName or 'endmodule' in endName:
                    if endName[0] == '!':
                        endName = '! Removed comment line'

            elif len(newSplits) > 1:
                endName = newSplits[0] + ' ' + newSplits[1]
                if endNameTarget in endName or endNameTarget1 in endName or 'end module' in endName or 'endmodule' in endName:
                    if endName[0] == '!':
                        endName = '! Removed comment line'

            else:
                endName = ''

        else:
            endName = ''

    if new:
        if 'end module' not in ' '.join(new.strip().lower().split()) or 'endmodule' not in ' '.join(new.strip().lower().split()):
            output.write('! Subprogram not used ' + new)

        else:
            output.write(new)

    try:
        innerStarts.pop(0)

    except IndexError:
        pass

    return innerStarts

for file in glob.glob(covPath + '/' + '*_UNCOVERED_FILES_COVERAGE_*'):
    soup = bs(open(file), 'html.parser')
    for l in soup.find_all('a'):

        if 'title' in l.attrs:
            delName = l.text
            
            if delName.endswith('.F90') or delName.endswith('.f90'):

                print delName

                try:
                    os.remove(srcPath + '/' + delName)

                except OSError:
                    pass


unused = {}
for file in glob.glob(covPath + '/' + '*_UNCVRDFNCS_*'):
    soup = bs(open(file), 'html.parser')
    modName = os.path.basename(soup.title.string.split()[-1])

    if modName not in unused:
        lnums = []
        print file

        for l in soup.find_all('a'):

            ltitle = l.attrs['title']
            if ltitle.startswith('Line'):
                lnums.append((int(ltitle.split()[-1]), l.text.lower()))

        unused[modName] = sorted(lnums)

#something is wrong with subcol_utils processing
ordered = ['subcol_allocate_internal', 'subcol_utils_init', 'subcol_get_nsubcol', 'subcol_get_indcol', 'subcol_get_filter', \
            'subcol_get_weight', 'subcol_get_ncol', 'subcol_set_nsubcol', 'subcol_set_filter', 'subcol_set_weight', \
            'is_weight_set', 'is_filter_set', 'is_subcol_on', 'subcol_get_scheme', 'subcol_utils_readnl', \
            'subcol_field_copy_1dint', 'subcol_field_copy_2dint', 'subcol_field_copy_3dint', 'subcol_field_copy_4dint', \
            'subcol_field_copy_5dint', 'subcol_field_copy_1ddouble', 'subcol_field_copy_2ddouble', 'subcol_field_copy_3ddouble', \
            'subcol_field_copy_4ddouble', 'subcol_field_copy_5ddouble', 'subcol_field_copy_1dreal', 'subcol_field_copy_2dreal', \
            'subcol_field_copy_3dreal', 'subcol_field_copy_4dreal', 'subcol_field_copy_5dreal', 'subcol_state_field_copy_1ddouble', \
            'subcol_state_field_copy_2ddouble', 'subcol_state_field_copy_3ddouble', 'subcol_state_copy', 'subcol_state_hdrinit', \
            'subcol_field_avg_shr_1dint', 'subcol_field_avg_shr_2dint', 'subcol_field_avg_shr_1ddouble', \
            'subcol_field_avg_shr_2ddouble', 'subcol_field_avg_shr_1dreal', 'subcol_field_avg_shr_2dreal', \
            'subcol_field_get_firstsubcol_1dint', 'subcol_field_get_firstsubcol_2dint', 'subcol_field_get_firstsubcol_1ddouble', \
            'subcol_field_get_firstsubcol_2ddouble', 'subcol_field_get_firstsubcol_1dreal', 'subcol_field_get_firstsubcol_2dreal', \
            'subcol_pack_1d_int', 'subcol_pack_2d_int', 'subcol_pack_3d_int', 'subcol_pack_4d_int', 'subcol_pack_5d_int', \
            'subcol_pack_6d_int', 'subcol_pack_1d_double', 'subcol_pack_2d_double', 'subcol_pack_3d_double', \
            'subcol_pack_4d_double', 'subcol_pack_5d_double', 'subcol_pack_6d_double', 'subcol_pack_1d_real', \
            'subcol_pack_2d_real', 'subcol_pack_3d_real', 'subcol_pack_4d_real', 'subcol_pack_5d_real', \
            'subcol_pack_6d_real', 'subcol_unpack_1d_int', 'subcol_unpack_2d_int', 'subcol_unpack_3d_int', \
            'subcol_unpack_4d_int', 'subcol_unpack_5d_int', 'subcol_unpack_6d_int', 'subcol_unpack_1d_double', \
            'subcol_unpack_2d_double', 'subcol_unpack_3d_double', 'subcol_unpack_4d_double', 'subcol_unpack_5d_double', \
            'subcol_unpack_6d_double', 'subcol_unpack_1d_real', 'subcol_unpack_2d_real', 'subcol_unpack_3d_real', \
            'subcol_unpack_4d_real', 'subcol_unpack_5d_real', 'subcol_unpack_6d_real', 'subcol_avg_inter_int', \
            'subcol_avg_inter_double', 'subcol_avg_inter_real', 'subcol_avg_int', 'subcol_avg_double', 'subcol_avg_real', \
            'subcol_utils_init_restart', 'subcol_utils_write_restart', 'subcol_utils_read_restart']

subcolset = set()
neworder = []
if 'subcol_utils.F90.in' in unused:
    unused['subcol_utils.F90'] = unused['subcol_utils.F90.in']
    for uu in unused['subcol_utils.F90']:
        subcolset.add(uu[1].strip().lower().split('_mp_')[-1][:-1])

    for o in ordered:
        if o.strip().lower() in subcolset:
            neworder.append((0, o.strip().lower()))

    unused['subcol_utils.F90'] = neworder


offsets = defaultdict(set)
subprogTypes = set(['integer', 'real(r8)', 'real(r4)', 'logical', 'integer(esmf_kind_i8)',\
                    'real(sp)', 'real(dp)', 'real(fp)', 'type(vdiff_selector)', \
                    'integer(shr_kind_i8)', 'real(shr_kind_r8)', 'real(kind=r8)', \
                    'character(len=9)', 'character*3', 'integer(i4)', 'character(len=16)', \
                    'integer(i4)'])
modifiers = set(['elemental', 'pure', 'recursive'])
for file, spLines in unused.iteritems():
    try:

        if file.endswith('.F90') or file.endswith('.f90'):

            modname = file.split('.')[0].lower()

            walk = readF(srcPath + '/' + file)
            starts = copy.deepcopy(spLines)
            with open('tmp', 'w') as output:
                for l in walk:

                    if starts:
                        lineSpace = ' '.join(l.strip().lower().split())
                        splitLine = lineSpace.split()
                        name = starts[0][1]
                        #processedName = name.replace(modname + '_mp_', '')[:-1]
                        if file != 'subcol_utils.F90':
                            processedName = name.split('_mp_')[-1][:-1]

                        else:
                            processedName = name

                        if splitLine and splitLine[0][0] != '!':
                            subprogType = splitLine[0]
                            if subprogType == 'subroutine':

                                if splitLine[1].split('(')[0] == processedName:
                                    spType = 'subroutine'
                                    starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                else:
                                    output.write(l)

                            elif subprogType == 'function':
                                if splitLine[1].split('(')[0] == processedName:
                                    spType = 'function'
                                    starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                else:
                                    output.write(l)

                            elif subprogType in subprogTypes or subprogType in modifiers or 'character(len=' in subprogType:
                                if splitLine[1] == 'subroutine':

                                    if splitLine[2].split('(')[0] == processedName:
                                        spType = 'subroutine'
                                        starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                    else:
                                        output.write(l)

                                elif splitLine[1] == 'function':
                                    if splitLine[2].split('(')[0] == processedName:
                                        spType = 'function'
                                        starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                    else:
                                        output.write(l)

                                elif splitLine[1] in modifiers:
                                    if splitLine[2] == 'subroutine':

                                        if splitLine[3].split('(')[0] == processedName:
                                            spType = 'subroutine'
                                            starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                        else:
                                            output.write(l)

                                    elif splitLine[2] == 'function':
                                        if splitLine[3].split('(')[0] == processedName:
                                            spType = 'function'
                                            starts = commentSubprog(walk, l, output, spType, processedName, starts)

                                        else:
                                            output.write(l)                                    

                                else:
                                    output.write(l)

                            else:
                                output.write(l)

                        else:
                            output.write(l)

                    else:
                        output.write(l)
                        
            os.rename('tmp', srcPath + '/' + file)

    except IOError:
        #print file
        pass
