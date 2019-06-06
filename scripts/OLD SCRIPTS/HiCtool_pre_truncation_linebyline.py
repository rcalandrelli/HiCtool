# Perform pre-truncation on reads that contain potential ligation junctions. To be executed before the mapping step.

# Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]
# Options:
#  -h, --help               show this help message and exit
#  -i INPUTFILES            Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].
#  -e RESTRICTION_ENZYMES   Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses a combination of MboI and Hinfl.

# Output files:
#  Fastq files with pre-truncated reads.
#  Log files with pre-truncation information.
#  Plot of the distribution of truncated reads length.

from optparse import OptionParser
from collections import Counter
import os
import re
from math import ceil
import matplotlib.pyplot as plt
from matplotlib import rcParams
from time import gmtime, strftime

parameters = {'inputFiles': None,
              'restriction_enzymes': None}

#class pre_truncate:
#    def __init__(self, parameters):
#        self.pre_truncation(parameters)

def pre_truncation(parameters):

    re_info = dict() # rs: restriction site sequence; lj: ligation junction sequence
    re_info['HindIII']={'rs':'AAGCTT','lj':'AAGCTAGCTT'}
    re_info['MboI']={'rs':'GATC','lj':'GATCGATC'}
    re_info['DpnII']={'rs':'GATC','lj':'GATCGATC'}
    re_info['Sau3AI']={'rs':'GATC','lj':'GATCGATC'}
    re_info['BglII']={'rs':'AGATCT','lj':'AGATCGATCT'}
    re_info['NcoI']={'rs':'CCATGG','lj':'CCATGCATGG'}
    re_info['Hinfl']={'rs':'GANTC','lj':'GA[ACGT]TA[ACGT]TC'}
    
    inputFiles = map(str, parameters['inputFiles'].strip('[]').split(','))
    restriction_enzymes = map(str, parameters['restriction_enzymes'].strip('[]').split(','))
    
    # Check that all the restriction enzymes are available
    for i in restriction_enzymes:
        if i not in re_info.keys():
            print "ERROR! " + i + " is not among the available restriction enzymes (HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl)! Check the spelling or contact Riccardo Calandrelli at <rcalandrelli@eng.ucsd.edu>."
            return
    
    for input_fastq in inputFiles:
        filename = os.path.splitext(input_fastq)[0]
        if len(restriction_enzymes) == 1:
            print "Pre-truncation of " + filename + ".fastq (restriction enzyme: " + restriction_enzymes[0] + ")..."
        else:
            print "Pre-truncation of " + filename + ".fastq (restriction enzymes: " + ", ".join(restriction_enzymes) + ")..."
        print "Start: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())
        
        count = 0 # to count reads with potential ligation junction
        lengths = [] # to save the length of the pieces after truncation
        n_reads = sum(1 for x in open(input_fastq))/4
        percents = {}
        for n in range(0,101,10)[1:]:
            percents[str(n) + '%'] = int(n_reads * n * 0.01)
        line_index = -1 # index to select the lines of the fastq file 4 by 4
        
        with open (input_fastq, 'r') as infile:
            for line in infile: # iteration over lines of the fastq file
                count += 1
                line_index += 1
                for key, value in percents.iteritems():
                    if line_index == value:
                        print key + ' completed.'
                if line_index == 1: # to save the read length to be used below for log files
                    read_length = len(line[:-1])
                if line_index % 4 == 0:
                    read_item = [] # list to save 4 lines of the fastq file per time re-initialized every time
                    read_item.append(line)
                elif line_index % 4 == 1 or line_index % 4 == 2:
                    read_item.append(line)
                elif line_index % 4 == 3:
                    read_item.append(line)
                    # Work on this block (4 lines) and search for the ligation junction
                    for j in restriction_enzymes:
                        match = re.search(re_info[j]['lj'],read_item[1])
                        if match: # ligation junction present in the read sequence
                            count += 1
                            sequence = read_item[1][0:-1] # remove the \n char at the end of the sequence
                            pieces = []
                            pieces.append(sequence[:match.start()])
                            pieces.append(sequence[match.end():])
                            max_length = max(len(x) for x in pieces)
                            if j != 'Hinfl':
                                re_seq = re_info[j]['rs']
                            else:
                                re_seq = match.group()[:4] + 'C'
                            lengths.append(max_length + len(re_seq))
                            piece = [x for x in pieces if len(x) == max_length][0]
                            start_index = re.search(piece,sequence).start()
                            match_index = match.start()
                            if start_index < match_index:
                                piece = piece + re_seq
                                read_item[1] = piece + '\n'
                                read_item[3] = read_item[3][:-1][start_index:start_index+len(piece)] + '\n'
                            elif start_index > match_index:
                                piece = re_seq + piece
                                read_item[1] = piece + '\n'
                                read_item[3] = read_item[3][:-1][start_index-len(re_seq):start_index-len(re_seq)+len(piece)] + '\n'
                    # Save the 4 lines to the output file
                    if line_index == 3:
                        with open (filename + '.trunc.fastq' ,'w') as fout:
                            for i in xrange(len(read_item)):
                                fout.write(read_item[i])
                    else:
                        with open (filename + '.trunc.fastq' ,'a') as fout:
                            for i in xrange(len(read_item)):
                                fout.write(read_item[i])
                        
        print '100% completed.'
            
        perc_reads = ceil(float(count)/float(n_reads)*10000)/100.0
        result = str(n_reads) + ' reads (length = ' + str(read_length) + ' bp); of these:\n  ' + str(count) + ' (' + str(perc_reads) + '%) contained a potential ligation junction and have been truncated.'
        print result

        with open ('pre_truncation_log.txt', 'a') as fout:
            if len(restriction_enzymes) == 1:
                fout.write(input_fastq + ', ' + restriction_enzymes[0] + '\n' + result + '\n\n')
            else:
                fout.write(input_fastq + ', ' + ", ".join(restriction_enzymes) + '\n' + result + '\n\n')
                
        rcParams.update({'figure.autolayout': True})
        cnt = Counter(lengths)
        my_lengths = [k for k, v in cnt.iteritems() if v >= 1]
        plt.close("all")
        plt.hist(lengths, bins=len(my_lengths))
        plt.title('Truncated reads (' + input_fastq + ')', fontsize=14)
        plt.xlabel('Read length', fontsize=12)
        plt.ylabel('Number of reads', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.savefig(input_fastq + '_truncated_reads.pdf', format = 'pdf')
        
        print "End: " + strftime("%Y-%m-%d %H:%M:%S", gmtime())


if __name__ == '__main__':
    
    usage = 'Usage: python2.7 HiCtool_pre_truncation.py [-h] [options]'
    parser = OptionParser(usage = 'python2.7 %prog -i inputFiles -e restriction_enzymes')
    parser.add_option('-i', dest='inputFiles', type='string', help='Input fastq file or input fastq files passed like a Python list. Example: [file1.fastq,file2.fastq].')
    parser.add_option('-e', dest='restriction_enzymes', type='string', help='Restriction enzyme (or restriction enzymes passed like a Python list: [enzyme1,enzyme2]). Choose between: HindIII, MboI, DpnII, Sau3AI, BglII, NcoI, Hinfl. Arima kit uses both MboI and Hinfl.')
    (options, args) = parser.parse_args( )

    if options.inputFiles == None:
        parser.error('-h for help or provide the input fastq file(s)!')
    else:
        pass
    if options.restriction_enzymes == None:
        parser.error('-h for help or provide the restriction enzyme(s)!')
    else:
        pass

    parameters['inputFiles'] = options.inputFiles
    parameters['restriction_enzymes'] = options.restriction_enzymes

    pre_truncation(parameters)