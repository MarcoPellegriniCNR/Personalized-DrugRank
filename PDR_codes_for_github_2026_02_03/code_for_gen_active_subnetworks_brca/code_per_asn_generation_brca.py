# codice epr genrare actrive subnetworks (modualre) per Breast cancer
#


filename_cp_modules = 'Core&Peel.txt'
filename_up_deg = 'TCGA_BRCA_population_UP_gene_list.txt'
filename_down_deg = 'TCGA_BRCA_population_DOWN_gene_list.txt'
filename_up_edger = 'TCGA_BRCA_population_UP_gene_edgeR.txt'
filename_down_edger = 'TCGA_BRCA_population_DOWN_gene_edgeR.txt'


# lettura del file dei moduli prodotto da core&peel


def get_modules(in_filename) :

    lines = []
    lines_as_lists = []
    total_set = set()
    with open(in_filename, encoding="utf-8") as f:
        line = f.readline()
        while line != '' :            
            lines.append(line.strip())
            line_array = line.split(' ')
            lines_as_lists.append(line_array)
            total_set.update(line_array)
            if len(line_array) > 100 :
                print('Exisitng long module', len(line_array))
            line = f.readline()   
        f.close()
        
    print('Numero moduli', len(lines))
    print('Numero geni unici', len(total_set))
    ##print(lines[0])
    ##print(lines[0].split(' '))

        
    return lines_as_lists, total_set

#modules, gene_set = get_modules(filename_cp_modules)

# test the crc asn forr long modules
#get_modules('CRC_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt')
#ok ci sono anceh i moduli lunghi

##############################

def get_up_down_genes() :

    with open(filename_up_deg, encoding="utf-8") as f:
        up_deg = f.read()
        up_genes_list = up_deg.split('\n')
        f.close()

    print('Number of up deg genes', len(up_genes_list),up_genes_list[0], up_genes_list[-1], up_genes_list[-2] )    

    with open(filename_down_deg, encoding="utf-8") as f:
        down_deg = f.read()
        down_genes_list = down_deg.split('\n')
        f.close()

    print('Number of down  deg genes', len(down_genes_list), down_genes_list[0], down_genes_list[-1], down_genes_list[-2] )    

    return(up_genes_list, down_genes_list )

#get_up_down_genes()

########################################


def get_up_down_edgeR_data(filename, lable):

    lines = []
    lines_as_lists = []
    print('Opening :', filename)
    with open(filename, encoding="utf-8") as f:
        line1 = f.readline()
        line1_array = line1.strip().split('\t')
        print(line1_array)
        line =  f.readline()
        while line != '' :            
            lines.append(line.strip())
            line_array = line.split('\t')
            lines_as_lists.append(line_array)
            assert len(line_array) == len(line1_array)
            line = f.readline()

        f.close()
    print('Number of edger records' + lable, len(lines_as_lists))
    return lines_as_lists
          


#up_edger_data = get_up_down_edgeR_data(filename_up_edger, 'up')
#down_edger_data = get_up_down_edgeR_data(filename_down_edger, 'down')

###################################################################


from pygeneconverter import ensembl_to_hugo, hugo_to_ensembl

def pygeneconverter_test():
    ensembl_list = ['ENSG00000223972', 'ENSG00000237613', 'ENSG00000268020']
    hugo_list = ['DDX11L1', 'FAM138A', 'OR4G4P']

    ensembl_df = ensembl_to_hugo(ensembl_list)
    hugo_df = hugo_to_ensembl(hugo_list)

    print(list(ensembl_df['HGNC_ID']), len(ensembl_df), ensembl_list)
    print(list(hugo_df['ENSEMBL_ID']), len(hugo_df),  hugo_list)


#pygeneconverter_test()
    

###############
def venn(set1, set2) :

    print('len1 {}, len2 {}, intersection {}, diff1-2 {}, diff2-1 {}, union {}'.format(len(set1), len(set2), len(set1 & set2) ,len(set1 - set2), len(set2 - set1), len(set1 | set2)))


########################################

    
#####################################################################

def pygeneconverter_application():
    filtered_up_edger_data = [x for x in up_edger_data if float(x[0]) >= 1.0 ]
    print('Filtered up genes', len(filtered_up_edger_data), filtered_up_edger_data[:10] )

    filtered_down_edger_data = [x for x in down_edger_data if float(x[0]) <= -1.0 ]
    print('Filtered down genes', len(filtered_down_edger_data), filtered_down_edger_data[:10])


    ######################################################################
    # convert to hugo names

    genes_up_edger_data = [x[5].strip() for x in filtered_up_edger_data ]
    a1 = ensembl_to_hugo(genes_up_edger_data)
    b1 = a1['HGNC_ID']
    c1 = list(b1)
    print(len(a1), len(b1), len(c1))
    hugo_genes_up_edger_data = c1

    print('Filtered up genes after mapping', len(hugo_genes_up_edger_data), hugo_genes_up_edger_data.count(''), hugo_genes_up_edger_data.count(None))

    genes_down_edger_data = [x[5].strip() for x in filtered_down_edger_data ]
    a2 = ensembl_to_hugo(genes_down_edger_data)
    b2 = a2['HGNC_ID']
    c2 = list(b2)
    print(len(a2), len(b2), len(c2))
    hugo_genes_down_edger_data =c2
    print('Filtered down genes after mapping', len(hugo_genes_down_edger_data), hugo_genes_down_edger_data.count(''), hugo_genes_down_edger_data.count(None))

    venn(gene_set, set(hugo_genes_up_edger_data + hugo_genes_down_edger_data))

    return gene_set &  set(hugo_genes_up_edger_data + hugo_genes_down_edger_data)

###########################################
# usng the hypergeoemtric test to filter modules.

deg_genes = pygeneconverter_application()
#
from scipy.stats import hypergeom

def apply_hpergeom_filter() :
    M = 12919 # total number of genes
    n =  2757 # toal number of BRCA deg genes in the M genes

    active_modules = []
    for module in modules :
        N = len(module)
        k = len(set(module) & deg_genes)
        p = hypergeom.sf(k-1, M, n, N)
        
        if p <= 0.01 :
    #        print (M,n, N, k, p)
            active_modules.append(module)

    print('Number of active modules' , len(active_modules))

    with open('BRCA_Core&Peel-coexpr025marta-r1-nl10-d0.7-f1-j0.8-marklargest_MultipleModules.txt', 'w') as f :
        for module in active_modules :
            f.write(' '.join(module))
        f.close()

    
# apply_hpergeom_filter() 
