"""
Program to:
1) Generate the global matrix containing all the contact matrices (intra and inter) for all the chromosomes to be normalized using run_ic_mes.sh (single processor).
2) Extract a single contact matrix from a global matrix (observed or normalized).
3) Plot the global contact matrix or a single contact matrix.

To use this code, an HiC_project_object.hdf5 must be provided (see HiCtool_hifive.py)
"""

chromosomes = {'hg38':{'1':248956422,
                   '2':242193529,
                   '3':198295559,
                   '4':190214555,
                   '5':181538259,
                   '6':170805979,
                   '7':159345973,
                   '8':145138636,
                   '9':138394717,
                   '10':133797422,
                   '11':135086622,
                   '12':133275309,
                   '13':114364328,
                   '14':107043718,
                   '15':101991189,
                   '16':90338345,
                   '17':83257441,
                   '18':80373285,
                   '19':58617616,
                   '20':64444167,
                   '21':46709983,
                   '22':50818468,
                   'X':156040895,
                   'Y':57227415},
                   'mm10':{'1':195471971,
                   '2':182113224,
                   '3':160039680,
                   '4':156508116,
                   '5':151834684,
                   '6':149736546,
                   '7':145441459,
                   '8':129401213,
                   '9':124595110,
                   '10':130694993,
                   '11':122082543,
                   '12':120129022,
                   '13':120421639,
                   '14':124902244,
                   '15':104043685,
                   '16':98207768,
                   '17':94987271,
                   '18':90702639,
                   '19':61431566,
                   'X':171031299,
                   'Y':91744698}}  


def save_matrix_rectangular(a_matrix, output_file):
    """
    Format and save an inter-chromosomal contact matrix in a txt file. 
    1) Data are reshaped to form a vector.
    2) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    3) Data are saved to a txt file.
    Parameters:
    a_matrix (numpy matrix): input contact matrix to be saved
    output_file: output file name in txt format
    Output:
    txt file containing the formatted data
    """
    import numpy as np
    n_row = np.shape(a_matrix)[0]
    n_col = np.shape(a_matrix)[1]
    vect = np.reshape(a_matrix,[1,n_row*n_col]).tolist()[0]
    with open (output_file,'w') as fout:
        k = len(vect)
        i = 0
        count = 0
        flag = False # flag to set if the end of the vector has been reached
        while i < k and flag == False:
            if vect[i] == 0:
                count+=1
                if (i+count == k):
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    flag = True
                    break
                while vect[i+count] == 0 and flag == False:
                    count+=1
                    if (i+count == k):
                        w_out = str(0) + str(count)
                        fout.write('%s\n' %w_out)
                        flag = True
                        break
                if flag == False:
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    i+=count
                    count = 0
            else:
                fout.write('%s\n' %vect[i])
                i+=1 


def load_matrix_rectangular(input_file, n_row, n_col):
    """
    Load a formatted contact matrix from a txt file and parse it.
    Parameters:
        input_file: input file name in txt format (generated by the function 
        "save_matrix_rectangular")
        n_row (int): number of rows of the matrix
        n_col (int): number of columns of the matrix
    Returns: 
        output_matrix: array containing the parsed values stored 
        in the input txt file to build a contact matrix        
    """
    import numpy as np    
    
    print "Loading " + input_file + "..."
    with open (input_file,'r') as infile:
        matrix_vect = []        
        for i in infile:
            if i[0] == "0" and i[1] != ".":
                for k in xrange(int(i[1:-1])):
                    matrix_vect.append(0)
            else:
                j = i[:-1]            
                matrix_vect.append(float(j))
  
    output_matrix = np.reshape(np.array(matrix_vect),[n_row,n_col])
    print "Done!"
    return output_matrix
   

def save_matrix(a_matrix, output_file):
    """
    Format and save an intra-chromosomal contact matrix in a txt file. 
    1) The upper-triangular part of the matrix is selected (including the
    diagonal).
    2) Data are reshaped to form a vector.
    3) All the consecutive zeros are replaced with a "0" followed by the
    number of times zeros are repeated consecutively.
    4) Data are saved to a txt file.
    Parameters:
    a_matrix (numpy matrix): input contact matrix to be saved
    output_file: output file name in txt format
    Output:
    txt file containing the formatted data
    """
    import numpy as np
    n = len(a_matrix)
    iu = np.triu_indices(n)
    vect = a_matrix[iu].tolist()
    with open (output_file,'w') as fout:
        k = len(vect)
        i = 0
        count = 0
        flag = False # flag to set if the end of the vector has been reached
        while i < k and flag == False:
            if vect[i] == 0:
                count+=1
                if (i+count == k):
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    flag = True
                    break
                while vect[i+count] == 0 and flag == False:
                    count+=1
                    if (i+count == k):
                        w_out = str(0) + str(count)
                        fout.write('%s\n' %w_out)
                        flag = True
                        break
                if flag == False:
                    w_out = str(0) + str(count)
                    fout.write('%s\n' %w_out)
                    i+=count
                    count = 0
            else:
                fout.write('%s\n' %vect[i])
                i+=1     


def load_matrix(input_file):
    """
    Load a formatted contact matrix from a txt file and parse it.
    Parameters:
        input_file: input file name in txt format (generated by the function 
        "save_matrix")
    Returns: 
        output_matrix: array containing the parsed values stored 
        in the input txt file to build a contact matrix        
    """
    import numpy as np    
    
    print "Loading " + input_file + "..."
    with open (input_file,'r') as infile:
        matrix_vect = []        
        for i in infile:
            if i[0] == "0" and i[1] != ".":
                for k in xrange(int(i[1:-1])):
                    matrix_vect.append(0)
            else:
                j = i[:-1]            
                matrix_vect.append(float(j))
  
    k = len(matrix_vect)
    matrix_size = int((-1+np.sqrt(1+8*k))/2)
    
    iu = np.triu_indices(matrix_size)
    output_matrix_1 = np.zeros((matrix_size,matrix_size)) # upper triangular plus the diagonal
    output_matrix_1[iu] = matrix_vect
    
    diag_matrix = np.diag(np.diag(output_matrix_1)) # diagonal
    output_matrix_2 = np.transpose(output_matrix_1) # lower triangular plus the diagonal
    output_matrix = output_matrix_1 + output_matrix_2 - diag_matrix
    print "Done!"
    return output_matrix
    
    
def save_matrix_tab(input_matrix, output_filename):
    """
    Save a contact matrix in a txt file in a tab separated format. Columns are
    separated by tabs, rows are in different lines.
    Parameters:
    input_matrix (numpy matrix): input contact matrix to be saved
    output_filename: output file name in txt format
    Output:
    txt file containing the tab separated data
    """
    with open (output_filename, 'w') as f:
            for i in xrange(len(input_matrix)):
                row = [str(j) for j in input_matrix[i]]
                if i != len(input_matrix) - 1:
                    f.write('\t'.join(row) + '\n')
                else:
                    f.write('\t'.join(row))
                    
                    
def load_matrix_tab(input_filename):
    """
    Load a contact matrix saved in a tab separated format using the function
    "save_matrix_tab".
    Parameters:
    input_filename (str): input contact matrix to be loaded
    Returns:
    output_matrix: array containing the parsed values stored 
    in the input tab separated txt file to build a contact matrix
    """
    import numpy as np
    
    print "Loading " + input_filename + "..."
    with open (input_filename, 'r') as infile:
        lines = infile.readlines()
        temp = []
        for line in lines:
            row = [float(i) for i in line.strip().split('\t')]
            temp.append(row)
            
        output_matrix = np.array(temp)
    print "Done!"
    return output_matrix


def load_topological_domains(input_file):
    """
    Function to load the topological domains coordinates from txt file.
    Parameters:
        input_file: input file name generated with "calculate_topological_domains" in txt format.
    Returns:
        List of lists with topological domain coordinates.
    """
    import csv
    print "Loading topological domain coordinates..."
    with open(input_file, 'r') as f:
        reader = csv.reader(f, dialect='excel', delimiter='\t')
        topological_domains = []
        for row in reader:
            row_int = [int(x) for x in row]
            topological_domains.append(row_int)
        print "Done!"
        return topological_domains

 
def generate_intrachromosomal_observed_data(a_chr,
                                            bin_size,
                                            input_file='HiC_project_object.hdf5',
                                            species='hg38',
                                            custom_species_sizes={},
                                            save_file=False):
    """
    Generate an observed intrachromosomal contact matrix from HiC_project_object.hdf5
    Parameters:
        a_chr (str): chromosome number (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        save_file (bool): if true, save the observed contact data.
    """
    import hifive
    
    chromosome = 'chr' + a_chr
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'mb_' + 'observed_fend'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome + '_' + bin_size_str + 'kb_' + 'observed_fend'    
    
    if species in chromosomes.keys():
        end_pos = (chromosomes[species][a_chr]/bin_size)*bin_size
    else:
        end_pos = (custom_species_sizes[a_chr]/bin_size)*bin_size
            
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.cis_heatmap(chrom=chromosome,
                                  start=0,
                                  stop=end_pos,
                                  binsize=bin_size,
                                  arraytype='full',
                                  datatype='raw')
    
    observed = heatmap_raw[:,:,0]
    
    if save_file == True:
        save_matrix(observed, output_filename + '.txt')  
    return observed


def generate_interchromosomal_observed_data(chr_row,
                                            chr_col,
                                            bin_size,
                                            input_file='HiC_project_object.hdf5',
                                            species='hg38',
                                            custom_species_sizes={},
                                            save_file=False):
    """
    Generate an observed interchromosomal contact matrix from HiC_project_object.hdf5
    Parameters:
        chr_row (str): chromosome number for the rows (example for chromosome 1: '1').
        chr_col (str): chromosome number for the columns (example for chromosome 1: '1').
        bin_size (int): bin size in bp of the contact matrix.
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        save_file (bool): if true, save the observed contact data.
    """
    import hifive
    
    chromosome_row = 'chr' + chr_row
    chromosome_col = 'chr' + chr_col
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        output_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + 'mb_'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        output_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + 'kb_'    
    
    if species in chromosomes.keys():
        end_pos_row = (chromosomes[species][chr_row]/bin_size)*bin_size
        end_pos_col = (chromosomes[species][chr_col]/bin_size)*bin_size
    else:
        end_pos_row = (custom_species_sizes[chr_row]/bin_size)*bin_size
        end_pos_col = (custom_species_sizes[chr_col]/bin_size)*bin_size
            
    hic = hifive.HiC(input_file)
    heatmap_raw = hic.trans_heatmap(chromosome_row, chromosome_col, 
                                    start1=0, stop1=end_pos_row, 
                                    start2=0, stop2=end_pos_col,
                                    binsize=bin_size, 
                                    datatype='raw')
    
    observed = heatmap_raw[:,:,0]
    row = observed.shape[0]
    col = observed.shape[1]
    
    if save_file == True:
        row_str = str(row)
        col_str = str(col)
        output_filename = output_filename + row_str + 'x' + col_str + '_'
        save_matrix_rectangular(observed, output_filename + 'observed.txt')
    return observed


def compute_matrix_data_full_observed(input_file='HiC_project_object.hdf5',
                                      bin_size=1000000,
                                      species='hg38',
                                      custom_species_sizes={},
                                      sexual_chromosomes=[],
                                      save_each_matrix=False,
                                      save_tab=True):
    """
    Generate the full observed contact matrix for all the chromosomes from HiC_project_object.hdf5.
    The matrix will contain all the chromosomes in the rows and in the columns, each entry i,j is the contact
    matrix associated to the chromosomes in row i and column j (either intra or interchromosomal).
    Parameters:
        input_file (str): object containing learned correction parameters in .hdf5 format obtained with
        HiCtool_hifive.py (default: 'HiC_project_object.hdf5').
        bin_size (int): bin size in bp of the contact matrix.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        save_each_matrix (bool): if true, save each single contact matrix in formatted txt file.
        save_tab (bool): if true, save the full observed matrix in tab separated format. This is
        needed to proceed with normalization using HiCorrector.
    Returns:
        Global matrix in numpy array format.
    """
    import numpy as np
    #sep_tick = 1 # width of the grid lines in pts
    #global matrix_global
    #global matrix_global_plot
    
    if species in chromosomes.keys():
        chromosomes_list = [str(i) for i in range(len(chromosomes[species]) - 1)[1:]] + ['X', 'Y']
    else:
        if len(sexual_chromosomes) > 0:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) - len(sexual_chromosomes) + 1)[1:]]
            chromosomes_list += sexual_chromosomes
        else:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) + 1)[1:]]
    
    for chr_row in chromosomes_list:
        chromosome_row = 'chr' + chr_row
        if chromosome_row == 'chr1':
            intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
            matrix_full_line = intra
            #matrix_full_line_plot = intra
        
        for chr_col in chromosomes_list:
            chromosome_col = 'chr' + chr_col
            
            if chromosome_row == chromosome_col and chromosome_row == 'chr1':
                continue
        
            elif chromosome_row == chromosome_col and chromosome_row != 'chr1':
                intra = generate_intrachromosomal_observed_data(chr_row,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
                n_row = np.shape(intra)[0]
                #sep_col = np.zeros((n_row,sep_tick))-1
                matrix_full_line = np.concatenate((matrix_full_line,intra),axis=1)
                #matrix_full_line_plot = np.concatenate((matrix_full_line_plot,sep_col,intra),axis=1)
                    
            else:
                row = (chromosomes[species][chr_row]/bin_size)*bin_size/bin_size
                col = (chromosomes[species][chr_col]/bin_size)*bin_size/bin_size
                row_str = str(row)
                col_str = str(col)
    
                matrix_data_full = generate_interchromosomal_observed_data(chr_row,chr_col,bin_size,input_file,species,custom_species_sizes,save_each_matrix)
                n_row = np.shape(matrix_data_full)[0]
                #sep_col = np.zeros((n_row,sep_tick))-1
                
                if 'matrix_full_line' in locals():
                    matrix_full_line = np.concatenate((matrix_full_line,matrix_data_full),axis=1)
                    #matrix_full_line_plot = np.concatenate((matrix_full_line_plot,sep_col,matrix_data_full),axis=1)
                else:
                    matrix_full_line = matrix_data_full
                    #matrix_full_line_plot = matrix_data_full
                   
        #n_col = np.shape(matrix_full_line_plot)[1]
        #sep_row = np.zeros((sep_tick,n_col))-1
        
        if chromosome_row == 'chr1':
            matrix_global = matrix_full_line
            #matrix_global_plot = matrix_full_line_plot
        else:
            matrix_global = np.concatenate((matrix_global,matrix_full_line))
            #matrix_global_plot = np.concatenate((matrix_global_plot,sep_row,matrix_full_line_plot))
            
        del(matrix_full_line)
        #del(matrix_full_line_plot)
    
    if bin_size >= 1000000:
        bin_size_str = str(bin_size/1000000)
        my_filename = 'HiCtool_' + bin_size_str + 'mb_'
    elif bin_size < 1000000:
        bin_size_str = str(bin_size/1000)
        my_filename = 'HiCtool_' + bin_size_str + 'kb_'
    
    save_matrix(matrix_global, my_filename + 'matrix_global_observed.txt') # without grid
    #save_matrix(matrix_global_plot, my_filename + 'matrix_global_observed_to_plot.txt') # with grid
    
    # Save the full observed matrix into a tab separated format
    if save_tab == True:
        save_matrix_tab(matrix_global, my_filename + 'matrix_global_observed_tab.txt')
                
    with open ('info.txt', 'w') as f:
        f.write('Rows: ' + str(len(matrix_global)) + '\n')
        f.write('Rowsum (average matrix * rows): ' + str(int(np.mean(matrix_global) * len(matrix_global))))
    
    return matrix_global
                
                
def extract_single_map(input_global_matrix,
                       tab_sep,
                       chr_row,
                       chr_col,
                       species='hg38',
                       bin_size=1000000,
                       data_type='observed',
                       custom_species_sizes={},
                       sexual_chromosomes=[],
                       save_output=True,
                       save_tab=False):
    """
    Extract a single contact matrix for a pair of chromosomes from the full matrix.
    Parameters:
        input_global_matrix (str): full contact matrix passed as a string of the filename saved to file.
        tab_sep (bool): if "input_global_matrix" is passed with a filename, then this boolean 
        tells if the global matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intrachromosomal map is extracted.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting. "observed" or "normalized".
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        save_output (bool): if true, save the contact matrix in formatted txt file.
        save_tab (bool): if true, save the contact matrix in tab separated format.
    Return:
        Contact matrix in numpy array format.
    """            
    if species in chromosomes.keys():
        chromosomes_list = [str(i) for i in range(len(chromosomes[species]) - 1)[1:]] + ['X', 'Y']
        chr_dim = []
        for i in chromosomes_list:
            chr_dim.append(chromosomes[species][i]/bin_size) 
        d_chr_dim = {}
        for i in chromosomes_list:
            d_chr_dim[i] = chromosomes[species][i]/bin_size
    else:
        if len(sexual_chromosomes) > 0:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) - len(sexual_chromosomes) + 1)[1:]]
            chromosomes_list += sexual_chromosomes
        else:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) + 1)[1:]]
        chr_dim = []
        for i in chromosomes_list:
            chr_dim.append(custom_species_sizes[i]/bin_size)                
        d_chr_dim = {}
        for i in chromosomes_list:
            d_chr_dim[i] = custom_species_sizes[i]/bin_size
    
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    # Loading global matrix
    if tab_sep == False:
        full_matrix = load_matrix(input_global_matrix)
    else:
        full_matrix = load_matrix_tab(input_global_matrix)
    
    if chr_row == '1':
        row_start = 0
    else:
        row_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_row)-1]]
    row_end = row_start + d_chr_dim[chr_row]
    
    if chr_col == '1':
        col_start = 0
    else:
        col_start = d_chr_dim_inc[chromosomes_list[chromosomes_list.index(chr_col)-1]]
    col_end = col_start + d_chr_dim[chr_col]
    
    output_matrix = full_matrix[row_start:row_end,col_start:col_end]
    
    if chr_row == chr_col:
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'mb_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_' + bin_size_str + 'kb_' + data_type + '.txt'
        if save_output == True:
            save_matrix(output_matrix, my_filename)
    else:
        dim_row = str(d_chr_dim[chr_row])
        dim_col = str(d_chr_dim[chr_col])
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'mb_' + dim_row + 'x' + dim_col + '_' + data_type + '.txt'
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000)
            my_filename = 'HiCtool_' 'chr' + chr_row + '_chr' + chr_col + '_' + bin_size_str + 'kb_' + dim_row + 'x' + dim_col + '_' + data_type + '.txt'
        if save_output == True:
            save_matrix_rectangular(output_matrix, my_filename)
    
    if save_tab == True:
        save_matrix_tab(output_matrix, my_filename.split('.')[0] + '_tab.txt')
    
    return output_matrix
    

def plot_map(input_global_matrix,
             tab_sep=False,
             chr_row='',
             chr_col='',
             bin_size=1000000,
             data_type='observed',
             species='hg38',
             custom_species_sizes={},
             sexual_chromosomes=[],
             my_colormap=['white', 'red'],
             cutoff_type='perc',
             cutoff=99,
             max_color='#460000',
             my_dpi=2000,
             topological_domains='',
             domain_color='#0000ff',
             plot_histogram=False):
    """
    Plot a contact map, either global or single map. To plot the global matrix leave "chr_row" and
    "chr_col" as empty strings.
    Parameters:
        input_global_matrix (str): full contact matrix passed as a string of the filename saved to file.
        tab_sep (bool): if "input_global_matrix" is passed with a filename, then this boolean 
        tells if the global matrix was saved in tab separated format (True) or not (False).
        chr_row (str): chromosome in the rows of the output contact matrix.
        chr_col (str): chromosome in the columns of the output contact matrix. If chr_col is 
        equal to chr_row then the intrachromosomal map is extracted.
        species (str): 'hg38' or 'mm10' or any other species label in string format.
        bin_size (int): bin size in bp of the contact matrix.
        data_type (str): which kind of data type you are extracting. "observed" or "normalized".
        custom_species_sizes (dict): dictionary containing the sizes of the chromosomes
        of your custom species. The keys of the dictionary are chromosomes in string
        format (example for chromosome 1: '1'), the values are chromosome lengths as int.
        sexual_chromosomes (list): list of the sexual chromosomes (if present) in your
        custom species (example for chromosome X: 'X').
        my_colormap (str | list): colormap to be used to plot the data. 1) Use a string if you choose among any colorbar here 
        https://matplotlib.org/examples/color/colormaps_reference.html 2) Use a list of strings with colors if you want
        a custom colorbar. Example: ['white', 'red', 'black']. Colors can be specified also in this format: '#000000'.
        cutoff_type (str): to select a type of cutoff ('percentile' or 'contact_number') or plot the full range of the data (set the 
        parameter as 'none').
        cutoff (int): percentile to set a maximum cutoff on the number of contacts for the colorbar.
        max_color (str): to set the color of the bins with contact counts over "cutoff".
        my_dpi (int): resolution of the contact map in dpi.
        topological_domains (str): topological domains txt file to visualize domains on the heatmap. If empty string, no topological domains.
        domain_color (str): to set the color for topological domains on the heatmap.
        plot_histogram (bool): if True, plot the contact data distribution.
    """         
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np   
    
    if species in chromosomes.keys():
        chromosomes_list = [str(i) for i in range(len(chromosomes[species]) - 1)[1:]] + ['X', 'Y']
        chr_dim = []
        for i in chromosomes_list:
            chr_dim.append(chromosomes[species][i]/bin_size) 
        d_chr_dim = {}
        for i in chromosomes_list:
            d_chr_dim[i] = chromosomes[species][i]/bin_size
    else:
        if len(sexual_chromosomes) > 0:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) - len(sexual_chromosomes) + 1)[1:]]
            chromosomes_list += sexual_chromosomes
        else:
            chromosomes_list = [str(i) for i in range(len(custom_species_sizes) + 1)[1:]]
        chr_dim = []
        for i in chromosomes_list:
            chr_dim.append(custom_species_sizes[i]/bin_size)                
        d_chr_dim = {}
        for i in chromosomes_list:
            d_chr_dim[i] = custom_species_sizes[i]/bin_size
    
    d_chr_dim_inc = {}
    k=1
    for i in chromosomes_list:
        d_chr_dim_inc[i] = sum(chr_dim[:k])
        k+=1
    
    d_chr_label_pos = {} # label position for chromosomes in the global matrix
    k = 0 # to consider the pixel occupied by the grid added after
    for i in chromosomes_list:
        d_chr_label_pos[i] = d_chr_dim_inc[i] - d_chr_dim[i]/2 + k
        k+=1
    
    label_pos = []
    label_name = []
    for i in chromosomes_list:
        label_pos.append(d_chr_label_pos[i])
        label_name.append('chr' + i)
    label_pos = np.array(label_pos)
    label_name = tuple(label_name)
    
    if chr_row == '' and chr_col == '':
        
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'HiCtool_' + bin_size_str + '_' + data_type
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'HiCtool_' + bin_size_str + '_' + data_type     
        
        if tab_sep == False:
            matrix_data_full = load_matrix(input_global_matrix)
        else:
            matrix_data_full = load_matrix_tab(input_global_matrix)
        
        print "Plotting..."
        # Adding grid to separate chromosomes
        k=0
        for i in chromosomes_list[:-1]:
            matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=1)
            matrix_data_full = np.insert(matrix_data_full, d_chr_dim_inc[i]+k, -1, axis=0)
            k += 1
        
        row = np.shape(matrix_data_full)[0]
        col = np.shape(matrix_data_full)[1]
        
        output_vect = np.reshape(matrix_data_full,row*col,1)
        negative_indexes = np.where(output_vect==-1)
        output_vect[negative_indexes] = 0
        non_zero = np.nonzero(output_vect)
        
        if isinstance(my_colormap, list):
            my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
        elif isinstance(my_colormap, str):
            my_cmap = my_colormap
        
        plt.close("all")
        
        if cutoff_type == 'perc':
            perc = np.percentile(output_vect[non_zero[0]],cutoff)
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=14)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == 'contact':
            perc = cutoff 
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc , vmin=0)
            plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=14)
            cbar = plt.colorbar(extend='max')
            cbar.cmap.set_over(max_color)
        elif cutoff_type == 'none':
            plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
            plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=14)
            cbar = plt.colorbar()
            
        cbar.cmap.set_under('black')   
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.xticks(label_pos, label_name, rotation='vertical', fontsize = 6)
        plt.yticks(label_pos, label_name, fontsize = 6)
        plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        print "Done!"
    
    else:
        
        matrix_data_full = extract_single_map(input_global_matrix,tab_sep,chr_row,chr_col,species,bin_size,data_type,custom_species_sizes,sexual_chromosomes,False,False)
        print "Plotting..."
        row = np.shape(matrix_data_full)[0]
        col = np.shape(matrix_data_full)[1]
        
        row_str = str(row)
        col_str = str(col)
        
        chromosome_row = 'chr' + chr_row
        chromosome_col = 'chr' + chr_col        
        
        if bin_size >= 1000000:
            bin_size_str = str(bin_size/1000000) + 'mb'
            my_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + row_str + 'x' + col_str + '_' + data_type
        elif bin_size < 1000000:
            bin_size_str = str(bin_size/1000) + 'kb'
            my_filename = 'HiCtool_' + chromosome_row + '_' + chromosome_col + '_' + bin_size_str + '_' + row_str + 'x' + col_str + '_' + data_type
        
            # Update matrix values to plot topological domains
        if topological_domains != '':
            if bin_size != 40000:
                print "WARNING! To plot topological domains the bin size should be 40000"
                return
            if chr_row != chr_col:
                print "WARNING! To plot topological domains the matrix should be intrachromosomal"
                return
            domains = load_topological_domains(topological_domains)
            diag_index = np.diag_indices(len(matrix_data_full))
            for domain in domains:
                temp_start = domain[0]/40000
                temp_end = domain[1]/40000
                matrix_data_full[temp_start,temp_start:temp_end] = -1
                matrix_data_full[temp_start:temp_end,temp_end-1] = -1
                matrix_data_full[(diag_index[0][temp_start:temp_end],diag_index[1][temp_start:temp_end])] = -1
        
        output_vect = np.reshape(matrix_data_full,row*col,1)
        non_zero = np.nonzero(output_vect)
        
        if isinstance(my_colormap, list):
            my_cmap = LinearSegmentedColormap.from_list('mycmap', my_colormap)
        elif isinstance(my_colormap, str):
            my_cmap = my_colormap
        
        def format_e(n):
            a = '%e' % n
            return a.split('e')[0].rstrip('0').rstrip('.') + 'e' + a.split('e')[1]        
        
        plt.close("all")
        plt.gcf().subplots_adjust(left=0.15)
        
        if cutoff_type == 'perc':
            perc = np.percentile(output_vect[non_zero[0]],cutoff)
            if topological_domains == '':
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
                cbar.cmap.set_under(domain_color)
        elif cutoff_type == 'contact':
            perc = cutoff 
            if topological_domains == '':
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmax=perc, vmin=0)
                cbar = plt.colorbar(extend='max')
                cbar.cmap.set_over(max_color)
                cbar.cmap.set_under(domain_color)
        elif cutoff_type == 'none':
            if topological_domains == '':
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest')
                cbar = plt.colorbar()
            else:
                plt.imshow(matrix_data_full, cmap=my_cmap, interpolation='nearest', vmin=0)
                cbar = plt.colorbar()
                cbar.cmap.set_under(domain_color)
        
        plt.title(data_type + ' contact map (' + bin_size_str + ')', fontsize=14)
        cbar.ax.set_ylabel(data_type + ' contact counts', rotation=270, labelpad=20)
        plt.ylabel(chromosome_row + ' coordinate (bp)', fontsize=12)
        plt.xlabel(chromosome_col + ' coordinate (bp)', fontsize=12)
        ticks_row = (np.arange(0, row, row/4) * bin_size)
        format_ticks_row = [format_e(i) for i in ticks_row.tolist()]
        ticks_col = (np.arange(0, col, col/4) * bin_size)
        format_ticks_col = [format_e(i) for i in ticks_col.tolist()]
        plt.yticks(np.arange(0, row, row/4), format_ticks_row)
        plt.xticks(np.arange(0, col, col/4), format_ticks_col)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.savefig(my_filename + '.pdf', format = 'pdf', dpi=my_dpi)
        
        # Plot of the histogram
        if plot_histogram:
            histogram = []
            if chr_row == chr_col:
                n = len(matrix_data_full)
                k = 1
                for i in xrange(n):
                    row = matrix_data_full[i][k:]
                    for j in row:
                        histogram.append(j)
                    k += 1
            else:
                histogram = matrix_data_full.reshape((1,row*col)).tolist()[0]
                
            plt.close("all")
            histogram_bins = int(pow(len(histogram),0.3))
            plt.hist(histogram, bins=histogram_bins)
            plt.title(data_type + ' contact counts distribution', fontsize=18)
            plt.xlabel(data_type + ' contact counts', fontsize=16)
            plt.ylabel('Number of bins', fontsize=16)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.tight_layout()
            plt.savefig(my_filename + '_histogram.pdf', format = 'pdf')
        
        print "Done!"
