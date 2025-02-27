##########################################################################################################
#  Module includes preprocessing and visualization functions for regulon analysis alongside pySCENIC.    #
#  Functions in the order of use in standard workflow. Dependencies included in environment .yaml        #
##########################################################################################################

# Requirements (Tentative)
import pyscenic 
import networkx as nx
import pandas as pd 
import numpy as np 
import scanpy as sc 
import seaborn as sns 
from matplotlib import pyplot as plt

################################
# pySCENIC Accessory Functions #
################################

def make_regulons(adjacencies_fname):
    '''
    Generates regulon object with pySCENIC load_motifs function, and a .csv generated from the first step of pySCENIC CLI 'grn'.
    '''
    from pyscenic.utils import load_motifs
    # Read in the adjacencies file created from the grn step. 
    adjacencies = pd.read_csv(ADJACENCIES_FNAME, index_col=False, sep=',')
    
    # Add motif metadata to the dataframe 
    df_motifs = load_motifs(REGULONS_FNAME)
    
    # Convert the dataframe to a sequence of regulons. 
    regulons = df2regulons(df_motifs)
    
    return regulons

def find_regulon(regulons_object, regulon_of_interest,verbose=True):
    """
    Scans a regulons object to return basic information of the regulon of interest.
    Stores the index of queried regulon for downstream purposes.
    """
    index = None  # Initialize index outside the loop
    for i, regulon in enumerate(regulons_object):  # Use enumerate for explicit indexing
        if hasattr(regulon, "name") and regulon.name == regulon_of_interest:  # Check for name attribute
            n_genes = len(regulons_object[i].gene2weight.keys())
            if verbose:
                print(f"{regulon_of_interest} found at index {i}, and has {n_genes} genes.")
            index = i
            break  # Exit loop if found

    if index is None:
        if verbose:
            print(f"Regulon '{regulon_of_interest}' not found in the regulons object.")
    return index 

def get_inflection_cutoff(regulon_object,plot:bool=False,point_selection:int=0):
    '''
    Calculates inflection point of target gene weight distributions of regulons within pySCENIC generated regulon object. 
    Point selection defaults to the first, to be most strict. 
    Prerequisite includes a regulon object derived from make_regulons()
    '''
    scores = []
    for i, regulon in enumerate(regulon_object):
        values = list(regulon_object[i].gene2weight.values())
        for value in values:
            scores.append(value)
            
    # Calculate inflection points based on target gene scores' distribution
    kde_object = stats.kde.gaussian_kde(scores)
    x_vals = np.linspace(min(scores), max(scores), 100)
    density_vals = kde_object(x_vals)
    first_derivative = np.gradient(density_vals)
    second_derivative = np.gradient(first_derivative)
    inflection_indices = np.where(np.diff(np.sign(second_derivative)))[0]
    inflection_points = x_vals[inflection_indices]
    print(f'inflection average method: {(inflection_points[0] + inflection_points[1] / 2)}')
    print(f'whole dictionary: {inflection_points}')
    inflection_threshold = inflection_points[point_selection] #first inflection point is default. 
    if plot:
        plt.plot(x_vals,density_vals)
        plt.axvline(x=inflection_threshold, color='red', linestyle='--', label='Inflection Point')
        plt.title('Density Plot of Gene Weights')
    
    return inflection_threshold

def filter_regulons(regulon_object,lower_cutoff:int=50,complexity_cutoff:int=10):
    '''''
    Generates target gene weight distribution inflection-based filter for a regulon object. Returns a filtered regulon object.

    lower_cutoff = minimum allowed number of target genes weight for the regulon to be retained pre filtering. 
    complexity_cutoff = minimum allowed number of target genes for a regulon to be retained post filtering.
    '''''
    cutoff_filtered_indices = []
    lower_count = 0
    inflection_cutoff = get_inflection_cutoff(regulon_object)

    for i, regulon in enumerate(regulon_object):
        # Store regulon names
        name = regulon.name
        names.append(name)

        # Deconstruct the gene2weight dictionary 
        weights = regulon.gene2weight
        if len(weights) < lower_cutoff:
            cutoff_filtered_indices.append(i)
            lower_count += 1 
    filtered_regulons = [value for index, value in enumerate(regulon_object) if index not in cutoff_filtered_indices]
    print(f'Filtered {lower_count} regulons with less than {lower_cutoff} target genes. {len(filtered_regulons)} regulons retained.')

    for i, regulon in enumerate(filtered_regulons):
        name = regulon.name
        weights = regulon.gene2weight
        weight_keys = weights.keys()
        weight_values = weights.values() 

        # Reconstruct weights 
        reconstructed_weights = dict(zip(weight_keys, weight_values)) # This is the mutable version of gene2weight{}
        initial_length = len(reconstructed_weights)
        
        # Filter the mutable dictionary 
        qualifying_genes = []
        for gene in reconstructed_weights:
            # Filter genes with score less than inflection cutoff
            if reconstructed_weights[gene] > inflection_cutoff:
                qualifying_genes.append(gene)

        filtered_targets_dict = {key: reconstructed_weights[key] for key in qualifying_genes if key in reconstructed_weights}
        print(f'Filtered {initial_length - len(filtered_targets_dict)} genes from {regulon.name} regulon: {len(filtered_targets_dict)} genes remain.')

        # Remove genes that have too few genes post-filtering
        if len(filtered_targets_dict) < complexity_cutoff:
            print(f'{regulon.name} regulon has less than {complexity_cutoff} genes post-filtering, it will be removed.')
        else:
            # Reconstruct regulons list with the filtered gene dictionary
            final_regulons.append(ctxcore.genesig.Regulon(name=name,gene2weight=filtered_targets_dict,gene2occurrence={},transcription_factor=name))
            
    return final_regulons

def regulons_to_csv(regulon_object, output_name):
    'Loads a dataframe of regulon data, then exports to .csv on output_name path'
    all_regulons = []
    for i, regulon in enumerate(regulon_object):
      # Iterate through gene-weight pairs
      for gene, weight in regulon.gene2weight.items():
        # Create a dictionary row with individual gene and weight
        row = {'Regulon':regulon.name ,'Target Genes': gene, 'Weights': weight}
        all_regulons.append(row)

    # Create the final DataFrame from the list
    df = pd.DataFrame(all_regulons)
    df.to_csv(f'{output_name}.csv', index=False)

def sort_auc(auc_df,adata):
    '''
    Pass a sorted adata object and apply the sorting to AUCell dataframes. (or any dataframe)
    '''
    cells = list(adata.obs.index)
    auc_df_sorted = auc_df.loc[cells]
    
    return auc_df_sorted

################################
#### Downstream Functions ######
################################

def color_mapping(source_data, palette:str='hls'):
    'creates color mapping based on some source data. Returns dictionary with components: row_colors, handles, labels.' 
    n_colors = len(set(source_data))
    #print(f'{n_colors} colors needed. For {len(source_data)} entities.')
    palette = sns.color_palette(palette, n_colors)

    # Map cell type to color in palette
    color_map = dict(zip(set(source_data), palette))
    row_colors = [color_map[entity] for entity in source_data]
    
    handles = [mpatches.Patch(color=color) for color in color_map.values()]
    labels = list(color_map.keys())
    
    color_info = {}
    color_info['row_colors'] = row_colors
    color_info['handles'] = handles 
    color_info['labels'] = labels
    #print(len(color_info['row_colors']))

    return color_info

def subset_data_pair(auc_mtx,adata,var:str,selection:str,preprocessed:bool=False):
    '''
    Create a touple of corresponding adata and AUCell dataframe.
    subset to contain cells meeting the same adata.obs condition(s).
    '''
    if preprocessed == False:
        adata_subset = adata[adata.obs[var] == selection]
        auc_subset = auc_mtx.loc[list(adata_subset.obs.index)]
        sorted_pair = [adata_subset,auc_subset]
    else:
        auc_subset = auc_mtx.loc[list(adata.obs.index)]
        sorted_pair = [adata,auc_subset]
    
    return sorted_pair 

def concat_objects(pair1,pair2,config):
    # [0] = adata [1] = aucell 
    '''
    After creating pairs for both genotypes, run this to concatenate them together for z-score comparisons.
    '''
    shared_ages = list(set(pair1[0].obs['age']).intersection(set(pair2[0].obs['age'])))
    a1 = pair1[0][pair1[0].obs['age'].isin(shared_ages)]
    a2 = pair2[0][pair2[0].obs['age'].isin(shared_ages)]
    ac = anndata.concat([a1,a2],join="outer")
    c1 = pair1[1][pair1[1].index.isin(ac.obs.index)]
    c2 = pair2[1][pair2[1].index.isin(ac.obs.index)] 
    cs = [c1,c2] 
    cc = pd.concat(cs) 

    cpair = subset_data_pair(cc,ac,var='age',selection='mOSN',preprocessed=True)
    cpair[0] = format_metadata(cpair[0],config,sort=True,concat=True)
    cpair[1] = sort_auc(cpair[1],cpair[0])
    
    return cpair

def plot_enrichment(auc_mtx,adata,color_basis:str,save:bool,path:str,concat:bool=False,z:bool=True,exp:bool=False,age_colors:bool=True,bin_by:str='age',ordered:bool=False):
    '''
    Create basic enrichment plot of regulons with color-coded cells based on a column in adata.obs. 
    '''
    if age_colors:
        ac = color_mapping(adata.obs['age'],palette='pastel')
        arc = ac['row_colors']
        ah = ac['handles']
        al = ac['labels']
        
    c = color_mapping(adata.obs[color_basis],palette='tab10')
    rc = c['row_colors']
    h = c['handles']
    l = c['labels']
    if concat:
        g = color_mapping(adata.obs['genotype'],palette='pastel')
        gr = g['row_colors']
        gh = g['handles']
        gl = g['labels']
        if exp:
            adata_df = adata.to_df()
            adata_df = adata_df[critical_list]
            #adata_df = adata_df.apply(zscore)
            if ordered:
                g = sns.clustermap(adata_df.T, figsize=(10,15), xticklabels=False, row_cluster = False,col_cluster = False,
                                   yticklabels=False,col_colors=[rc,gr],vmin=-2,vmax=2)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h+gh,
                           labels = l+gl,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 
            else:
                g = sns.clustermap(adata_df.T, figsize=(10,15), xticklabels=False, row_cluster = True,col_cluster = False,
                                   yticklabels=False,col_colors=[rc,gr],vmin=-2,vmax=2)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h+gh,
                           labels = l+gl,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 
        else:
            
            if z:
                #apply zscore
                za = auc_mtx.apply(zscore)
                if ordered:
                    g = sns.clustermap(za.T, figsize=(10,15), xticklabels=False, row_cluster = False,col_cluster = False,
                                       yticklabels=True,col_colors=[rc,gr],vmin=-2,vmax=2)
                    g.ax_heatmap.grid(False)
                    plt.legend(handles= h+gh,
                               labels = l+gl,
                               title = "Legend",
                               ncol = 2,
                               bbox_to_anchor = (1,1), 
                               bbox_transform = plt.gcf().transFigure,
                               loc = "upper right") 
                else:
                    g = sns.clustermap(za.T, figsize=(10,15), xticklabels=False, row_cluster = True,col_cluster = False,
                                       yticklabels=True,col_colors=[rc,gr],vmin=-2,vmax=2)
                    g.ax_heatmap.grid(False)
                    plt.legend(handles= h+gh,
                               labels = l+gl,
                               title = "Legend",
                               ncol = 2,
                               bbox_to_anchor = (1,1), 
                               bbox_transform = plt.gcf().transFigure,
                               loc = "upper right") 
            
    else:
        if z:
            za = auc_mtx.apply(zscore) 
            if ordered:
                g = sns.clustermap(za.T, figsize=(10,15), xticklabels=False, row_cluster = False,col_cluster = False,
                                       yticklabels=True,col_colors=[arc,rc],vmin=-2,vmax=2)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h+ah,
                           labels = l+al,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 
            else:
                g = sns.clustermap(za.T, figsize=(10,15), xticklabels=False, row_cluster = True,col_cluster = False,
                                       yticklabels=True,col_colors=[arc,rc],vmin=-2,vmax=2)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h+ah,
                           labels = l+al,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 
        else:
            if ordered:
                g = sns.clustermap(auc_mtx.T, figsize=(10,15), xticklabels=False, row_cluster = False,col_cluster = False,
                                   yticklabels=True,col_colors=rc)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h,
                           labels = l,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 
            else:
                g = sns.clustermap(auc_mtx.T, figsize=(10,15), xticklabels=False, row_cluster = True,col_cluster = False,
                                   yticklabels=True,col_colors=rc)
                g.ax_heatmap.grid(False)
                plt.legend(handles= h,
                           labels = l,
                           title = "Legend",
                           ncol = 2,
                           bbox_to_anchor = (1,1), 
                           bbox_transform = plt.gcf().transFigure,
                           loc = "upper right") 

    if save:
        plt.savefig(f'{path}.pdf', format="pdf", transparent = True)

def sigTest_regulons(integrated_pair:list,selection:str,selection_var:str,group_var:str='genotype'):
    '''
    Performs significance thresholding on regulons' enrichment, for concatenated pair of auc matrices and adata. 
    Plots filtered auc matrix using 'plot_enrichment()'
    '''
    # Initialize necessary objects and subsets
    ip = integrated_pair
    #ugv = unique group variables.
    ugv = list(ip[0].obs[group_var].unique())
    print(f'Will compare {selection} cells between: {ugv}') 
    # z-score transform the non-subset auc matrix. 
    zauc = ip[1].apply(zscore)
    # This subsets auc mtx by cells within the desires selection. Still contains both groups. 
    subset_auc = zauc[zauc.index.isin(ip[0].obs[ip[0].obs[selection_var] == selection].index)]
    g1 = subset_auc[subset_auc.index.isin(ip[0].obs[ip[0].obs[group_var] == 'Wt'].index)]
    g2 = subset_auc[subset_auc.index.isin(ip[0].obs[ip[0].obs[group_var] == 'Fzd1KO'].index)]

    print(f'wt {len(g1)} cells parsed ,ko {len(g2)} cells parsed')
    regulon_names = list(g1.columns)
    passing_regulons = []
    for regulon in regulon_names:
        # wilcoxon rank sum test between columns
        stat, pval = stats.ttest_ind(g1[regulon],g2[regulon])
        if pval < 0.01:
            passing_regulons.append(regulon)
    print(f'{len(passing_regulons)} out of {len(regulon_names)} are significantly different.')
    # Visualize significance-subset enrichment.
    plot_enrichment(subset_auc[passing_regulons],ip[0],color_basis='age',concat=True,save=False,z=False)

def compare_regulons(integrated_pair:list,regulon_object,regulon_pair:list,export:bool,path:str):
    '''
    Finds intersecting genes between two regulons. Generates dataframe with a column per regulon, and a column with intersection.
    Export on path for example: '/home/analysis/supplementary_tables/' 
    '''
    ro = regulon_object
    r1 = regulon_pair[0]+'(+)'
    r2 = regulon_pair[1]+'(+)'
    r1_genes = list(ro[find_regulon(ro,r1)].gene2weight.keys())
    r2_genes = list(ro[find_regulon(ro,r2)].gene2weight.keys())
    intersection = np.intersect1d(r1_genes,r2_genes)
    results =[r1_genes,r2_genes,intersection]
    results = pd.DataFrame(results)
    results = results.T
    if export:
        results.to_csv(f'{path}{regulon_pair[0]}_{regulon_pair[1]}_compared.csv')

def grab_targets(adjacency_path, regulon_object,critical_targets):
    '''
    Given target genes of interest or "critical_targets", parse a regulon object and grn adjacency file 
    to compose a table of adjacencies and regulons which have connectivity to them.
    '''
    crits = {}
    # Parse adjacency file from pySCENIC. Extract genes connected to desired list of genes.
    adjacencies = {}
    regulons = {}
    tfs = {}

    # Build the adjacencies dataframe
    for target in critical_targets:
        # Find adjacencies for target in grn file.
        adjacencies[target] = []
        with open(adjacency_path, 'r') as f:
            for line in f:
                line = str(line).strip().split(',')
                source = line[0]
                dest = line[1]
                weight = line[2] 
                if dest == target:
                    adjacencies[target].append(source)
            f.close()
    
        # Build regulon target dataframe
        regulon_names = []
        for regulon in regulon_object:
            regulon_names.append(regulon.name)

        for regulon in regulon_object:
            tfs[regulon.name] = regulon.name.replace('(+)','') # Extract just gene name.
            regulons[regulon.name] = list(regulon.gene2weight.keys()) # Extract all regulon targets
            if target in regulon.gene2weight.keys(): # If the critical target is in the regulon's targets
                if target == regulon.name.replace('(+)',''): # Don't retain regulons which target themselves only. 
                    pass
                else:
                    if regulon.name not in crits.keys(): # If it doesn't exist in the dictionary already
                        crits[regulon.name] = [] # Initialize the list, then append. 
                        crits[regulon.name].append(target)
                    else: # If it exists, just append the value list. 
                        crits[regulon.name].append(target)
        reg_df = {'Regulon': regulon_names}
        reg_df = pd.DataFrame(reg_df)
        reg_df['Targets'] =  reg_df['Regulon'].map(regulons)
        reg_df['critical_targets'] = reg_df['Regulon'].map(crits)
        reg_df = reg_df.replace(np.nan,None)
        reg_df['TF'] = reg_df['Regulon'].map(tfs) 
        adj_df = {'Target': list(adjacencies.keys()), 'Source': list(adjacencies.values())}
        adj_df = pd.DataFrame(adj_df)

                    
    return {'Adj':adj_df,'Reg':reg_df}               

def enrichment_correction(auc_mtx,adata,scaling_factor:int=.001,partition:str='celltype'):
    '''
    For cases where you bin an enrichment score mtx (auc_mtx) , apply normalization to a matrix based on the number of cells within a bin.
    This removes bias in representative enrichment scores based on the number of cells in the bin, and should reduce effect of outliers. 
    Scaling factor is user-set for interpretability, and is statistically arbitrary. Default of 1000 is used. 

    Expects a raw AUCell matrix.
    Returns a corrected matrix. 

    (Only situationally useful! Likely for testing.)
    '''
    for val in adata.obs[partition].unique():
        print(f'Normalizing {val} cells')
        # Make vector of cells that meet a particular criteria. 
        indexes = list(adata.obs[adata.obs[partition] == val].index)
        auc_mtx = auc_mtx.T # Makes cell IDs columns
        non_zeroes = len(auc_mtx[indexes]!=0)
        auc_mtx[indexes] = (auc_mtx[indexes] * non_zeroes) * scaling_factor
        auc_mtx = auc_mtx.T
    return auc_mtx 

def bin_enrichment(auc_mtx,adata,basis:str,celltype_subset:bool=True,scaling_factor:float=0.001,method:str='mean'):
    '''
    Bins auc mtx based on an adata.obs variable. Each bin will represent the median enrichment of cells within it. (WIP)
    Returns z-scored and non-zscored binned dataframes. 
    '''
    auc = auc_mtx.apply(zscore)
    #auc = auc.apply(zscore)
    #print(auc.columns)
    sub_cells = {}
    if basis == 'celltype':
        print(f'Grouping t-OSN cells with their respective parental clusters.')
        adata = adata.copy()
        adata.obs['celltype'] = adata.obs['celltype'].replace('tOSN-i','iOSN')
        adata.obs['celltype'] = adata.obs['celltype'].replace('tOSN-m','mOSN')
    for val in list(adata.obs[basis].unique()):
        cells  = list(adata.obs[adata.obs[basis] == val].index)
        sub_cells[val] = cells 
    binned_results = {}
    for subset in sub_cells.keys():
        binned_results[subset] = {}
        mask = sub_cells[subset]
        auc_sub = auc.T[mask]
        for regulon, cells in auc_sub.iterrows():
            if method == 'mean':
                binned_results[subset][regulon] = np.mean(cells)
            if method == 'sum':
                binned_results[subset][regulon] = np.sum(cells) 
            if method == 'median':
                binned_results[subset][regulon] = np.median(cells)                 
            #print(f'{subset}:{regulon} mean = {binned_results[subset][regulon]}')
    binned_results = pd.DataFrame.from_dict(binned_results)
    #binned_results = binned_results.rank(axis=1) #if erroneous  
        
    return binned_results

def get_cell_distributions(adata): 
    '''
    Accessory function only situationally useful. For this analysis, we had stages where transitional clusters of OSNs were examined.
    This function consolidates them and provides basic information for age/celltype distributions given adata object. 
    '''
    adata = adata.copy()
    adata.obs['celltype'] = adata.obs['celltype'].replace('tOSN-i','iOSN')
    adata.obs['celltype'] = adata.obs['celltype'].replace('tOSN-m','mOSN')
    distributions = {}
    totals = {}
    for age in adata.obs['age'].unique():
        distributions[age] = {} 
        subset = adata.obs[adata.obs['age'] == age]
        total = len(subset.index)
        totals[age] = total
        for celltype in subset['celltype'].unique():
            distributions[age][celltype] = len(subset[subset['celltype']==celltype].index) / total
    distributions = pd.DataFrame.from_dict(distributions)
    distributions.fillna(0, inplace=True)
    distributions = distributions.iloc[[1,0,2,3,4]]
    totals = pd.DataFrame(totals,index=[0])
    return distributions, totals

def celltype_genotype_comparison(path:str,cs_dict:dict,auc1,auc2,a1,a2):
    '''
    Pass a dictionary with cell-type as keys, and a list of specific regulons per cell type (determined via enrichment pattern, or RSS via pyscenic) as values (cs_dict),
    alongside auc matrices(auc1,auc2) for each genotype to be compared. This function calculates enrichment over time per cell-type 
    subset between conditions, and produces a heat-based subtraction showing how enrichment varies over time for a group of cell-type specific regulons. 
    '''
    # Consolidate tOSN celltypes 
    #a1 = a1[a1.obs['age'].isin(['P3','P5','P7','P10'])]
    #auc1 = auc1[auc1.index.isin(a1.obs.index)]
    subtractions = {}
    for celltype in cs_dict.keys(): 
        # For each celltype, generate two dataframes. (For subtraction) 
        for i, adata in enumerate([a1,a2]): 
            if i == 0:
                if celltype == 'mOSN':
                    sub = a1[a1.obs['celltype'].isin(['mOSN','tOSN-m'])]
                if celltype == 'iOSN':
                    sub = a1[a1.obs['celltype'].isin(['iOSN','tOSN-i'])]
                else:
                    sub = a1[a1.obs['celltype'] == celltype]
            if i == 1:
                if celltype == 'mOSN':
                    sub = a2[a2.obs['celltype'].isin(['mOSN','tOSN-m'])]
                if celltype == 'iOSN':
                    sub = a2[a2.obs['celltype'].isin(['iOSN','tOSN-i'])]
                else:
                    sub = a2[a2.obs['celltype'] == celltype]
                
            # 0 = wt , 1 = KO | subset auc to be the same cells as subset adata. 
            if i == 0:
                sub_auc = auc1[auc1.index.isin(sub.obs.index)]
                # Subset to celltype specific regulons
                sub_auc = sub_auc[cs_dict[celltype]]
                sub_binned = bin_enrichment(sub_auc,sub,basis='age')
                x = sns.clustermap(sub_binned,col_cluster=False,figsize= (5,5),cmap='Blues',vmin=-1,vmax=1)
                plt.suptitle(f'{celltype} Specific Enrichment Over Time: WT')
                x.ax_row_dendrogram.set_visible(False)
                plt.close()
                ordered_rows = x.dendrogram_row.reordered_ind
                original_binned = sub_binned.iloc[ordered_rows]
            if i == 1:
                sub_auc = auc2[auc2.index.isin(sub.obs.index)]
                sub_auc = sub_auc[cs_dict[celltype]]
                sub_binned = bin_enrichment(sub_auc,sub,basis='age')
                reordered_binned = sub_binned.iloc[ordered_rows]
                y = sns.clustermap(reordered_binned,col_cluster=False,row_cluster=False,figsize= (5,5),cmap='Blues',vmin=-1,vmax=1)
                plt.suptitle(f'{celltype} Specific Enrichment Over Time: KO')
                plt.close()
                subtraction = reordered_binned - original_binned
                subtractions[celltype] = subtraction
                z = sns.clustermap(subtraction,col_cluster=False,row_cluster=False,figsize=(5,5),cmap='vlag')
                plt.suptitle(f'{celltype} WT-KO Enrichment Over Time')    
                plt.savefig(f'{path}{celltype}_enr_subt_agebin.pdf', format="pdf", transparent = True)
    return subtractions            

def hash_variable(variable_name,adata):
    'generalization of metadata for mapping. appends salt prefix based on variable name to integer data, converts str data to hash.'
    variable_name = str(variable_name)
    salt = str(hash('1'))
    adata.obs['hash_key'] = [int(salt + str(hash(val)) if isinstance(val, int)
                             else hash(val)) for val in adata.obs[variable_name]]
    skeleton_key = dict(zip(adata.obs['hash_key'],adata.obs[variable_name]))
    return skeleton_key

def animate_umap(adata,frame_dir,color_by:str,calculate_umap:bool):
    '''
    Prepare 3D umap and frame directory for imagemagick and convert frames to .gif // Adapted from Sanbomics YouTube Channel.

    Limitations for Thursday:
        1. os call for running imagemagick doesn't work since it doesn't run from shell (apparently) 
        2. current color limitation is 20 for tab20 palette. Looking for alternative ways to add more colors and scale legend based on num_colors automatically.
        3. general efficiency and speed. (automatically remove the *.png from frame_directory after imagemagick converts to .gif. 
    '''
    if calculate_umap:
        print('Calculating 3D UMAP')
        sc.tl.umap(adata, n_components = 3)
        
    # Hash conversion of variable
    skeleton_key = hash_variable(color_by, adata)
    n_colors = len(adata.obs['hash_key'].unique())
    if n_colors <= 10:
        color = dict(zip(range(0,n_colors), plt.cm.tab10(range(0,n_colors))))
        columns = 2
    else:
        columns = 5
        color = dict(zip(range(0,n_colors), plt.cm.tab20(range(0,n_colors))))
    color = dict(zip(adata.obs['hash_key'].unique(),color.values()))
    umap = adata.obsm['X_umap']
    frame_dir = frame_dir

    # Redirect hash-dictionary to color arrays for easy legend creation
    legend_guide = {skeleton_key[key]: color[key] for key in color if key in skeleton_key}

    # Create legend with manual patches
    colors = adata.obs['hash_key'].map(color)
    legend_texts = list(legend_guide.keys())
    legend_colors = list(legend_guide.values())
    patches = [ mpatches.Patch(color = legend_colors[i], label ="{:s}".format(legend_texts[i]))for i in range(len(legend_texts)) ]

    # Plot frames for animation ( This section mostly comes from Sanbomics, I just add some QoL) 
    print('Generating 3D Frames')
    for i in tqdm(range(0,360,2), total = len(range(0,360,2))):
        fig = plt.figure(figsize = (10,10))
        plt.legend(handles= patches,
               title = "Legend",
               ncol = columns,
               bbox_to_anchor = (1,1), 
               bbox_transform = plt.gcf().transFigure,
               loc = "upper right") 
        ax = fig.add_subplot(projection = '3d')
        ax.scatter(umap[:,0],umap[:,1],umap[:,2], c = colors)

        x_center = (umap[:,0].max()+umap[:,0].min())/2
        y_center = (umap[:,1].max()+umap[:,1].min())/2
        z_center = (umap[:,2].max()+umap[:,2].min())/2

        ax.plot([x_center,x_center], [y_center,y_center], [umap[:,2].min() - 2, umap[:,2].max() + 2], c = 'k', lw = 5)
        ax.plot([x_center,x_center], [umap[:,1].min() -2, umap[:,1].max() + 2], [z_center,z_center], c = 'k', lw = 5)
        ax.plot([umap[:,0].min() -2, umap[:,0].max() +2],[y_center,y_center],[z_center,z_center],c = 'k', lw = 5)

        ax.view_init(20,i)
        ax.axis('off')

        plt.savefig(f'{frame_dir}/{i:003}.png',dpi = 100, facecolor = 'white')
        plt.close(fig)
        
    # Calls !imagemagick convert to create final .gif 
    print('All frames saved, running imagemagick on frames directory!')
    subprocess.run(f"convert -delay 5 {frame_dir}/*png {color_by}_umap.gif",stdout =subprocess.PIPE,shell = True)

def group_OSNs(adata,config):
    '''
    groups all OSN celltypes as one OSN cluster to compare to others.
    '''
    grouping = {} 
    for i, celltype in enumerate(config['Celltypes'].keys()):
        if celltype in ['tOSN-i','iOSN','tOSN-m','mOSN']:
            grouping[celltype] = 'OSN' 
        else:
            grouping[celltype] = celltype 
    adata.obs['OSN_class'] = adata.obs['celltype'].map(grouping)
    return adata 

def make_markers_unique(adata,var,premise:str,n_genes:int):
    '''
    Fixes sc.tl.rank_genes_groups to extract top n unique markers per group var. Sorts dataframe by a premise such as logfoldchanges or pval.
    '''
    ranking_dict = {}
    for group in adata.obs[var].unique(): 
        ranking_dict[group] = sc.get.rank_genes_groups_df(adata,group).sort_values(by=premise,ascending=False)

    top_genes = {}
    for df in ranking_dict:
        top_genes[df] = list(ranking_dict[df]['names'][:n_genes])
    if var == 'age':
        order = ['P3','P5','P10','P7','P14','P21','Adult','E13.5','E18.5','P0']
        top_genes = {key: top_genes[key] for key in order}
        fig, ax = plt.subplots(1,1, figsize = (15,10))
        sc.pl.rank_genes_groups_matrixplot(adata,
                                           values_to_plot='logfoldchanges',
                                           var_names=top_genes,
                                           ax = ax,
                                          vmin = -2,
                                          vmax = 2,
                                          dendrogram=False,
                                          categories_order = ['P3','P5','P10','P7','P14','P21','Adult','E13.5','E18.5','P0'],
                                          cmap='coolwarm',
                                          save=f'wt_{var}_binned_degs.pdf')
    if var == 'OSN_class':
        order = ['HBC','GBC','INP','OSN']
        top_genes = {key:top_genes[key] for key in order}
        fig, ax = plt.subplots(1,1, figsize = (15,10))
        sc.pl.rank_genes_groups_matrixplot(adata,
                                           values_to_plot='logfoldchanges',
                                           var_names=top_genes,
                                           ax = ax,
                                          vmin = -2,
                                          vmax = 2,
                                          dendrogram=False,
                                          categories_order = ['HBC','GBC','INP','OSN'],
                                          cmap='coolwarm',
                                          save=f'wt_{var}_binned_degs.pdf')
    return top_genes

def what_targets_gene(gene,target_df_dict,list_input:bool,verbose:bool=False):
    'Summarizes regulons and/or genes which target user input gene, given the output of multi_target_query()'
    results = {}
    if list_input:
        gene_list = gene
        for gene in gene_list:
            targets = []
            regulons = []
            results[gene] = {'Reg': [],'Adj': []}
            for row in target_df_dict['Reg'].iterrows():
                if gene in list(row[1]['Targets']):
                    regulons.append(row[1]['Regulon'])
            for row in target_df_dict['Adj'].iterrows():
                if gene == row[1]['Target']:
                    targets = row[1]['Source']
            if len(regulons) == 0:
                regulons = 'None'
            results[gene]['Reg'] = regulons
            results[gene]['Adj'] = targets 
            if verbose:
                if len(regulons) > 0:
                    print(f'{gene} is targetted by {regulons} regulons.')
                if len(targets) > 0:
                    print(f'{gene} is secondarily connected to: {targets}')
                if len(targets) + len(regulons) == 0:
                    print(f'Cannot derive what targets {gene} with current data.')
            
    else: 
        targets = []
        regulons = []
        results[gene] = {'Reg': [], 'Adj': []}
        for row in target_df_dict['Reg'].iterrows():
            if gene in list(row[1]['Targets']):
                regulons.append(row[1]['Regulon'])
        for row in target_df_dict['Adj'].iterrows():
            if gene == row[1]['Target']:
                targets = row[1]['Source']
        if len(regulons) == 0:
            regulons = 'None'
        results[gene]['Reg'] = regulons
        results[gene]['Adj'] = targets
        if verbose:
            if len(regulons) > 0:
                print(f'{gene} is targetted by {regulons} regulons.')
            if len(targets) > 0:
                print(f'{gene} is connected by the following genes as GRN adjacencies: {targets}')
            if len(targets) + len(regulons) == 0:
                print(f'Cannot derive what targets {gene} with current data.')
            
    #if list_input:  
        #results  = {'gene' : gene_list , 'Regulons': results[gene]['Reg'],'Adjacencies': results[gene]['Adj']}

    return results

def secondary_connection(gene, target_df_dict):
    'Looks at what_targets output and sees regulons targetting the adjacencies connected to query gene' 
    #print('-----Defining primary connections-----\n')
    root = what_targets_gene(gene,target_df_dict,list_input = False) 
    gene_list = root[gene]['Adj']
    #print('\n-----Defining secondary connections-----\n')
    secondary = what_targets_gene(gene_list,multi_targets,list_input = True)
    prim = pd.DataFrame.from_dict(root).T
    sec = pd.DataFrame.from_dict(secondary).T
    sec = sec.reset_index()
    sec = sec.rename(columns = {'index':f'{gene}_connections'})
    return [prim, sec]
    
def expression_from_network(G,scaled_ordered_result:list):
    reg = []
    genes = []
    # Prep lists to subset respective matrices.
    for node in G.nodes:
        if '(+)' in node:
            reg.append(node) 
        else:
            genes.append(node) 

    ssx = scaled_ordered_result[0][genes]
    print(ssx) 
    ssa = scaled_ordered_result[1][reg]

def unionize_network(root_gene:str,path:str,sec_df,prim_df):
    '''
    Draw network visualization of connections between primary and secondary targets to a gene of interest.

    sec_df and prim_df are output from secondary_connection()
    ''' 
    # Initialize graph
    G = nx.Graph()
    G.add_node(root_gene, color = 'blue') 
    for row, data in prim_df.iterrows(): # Attach root gene to establisged regulon connections. 
        regulon = data['Reg'] 
        if regulon == 'None':
            pass
        else: 
            for regulon in regulon:
                G.add_node(regulon, color = 'red')
                G.add_edge(regulon, root_gene) 
                
    for row, data in sec_df.iterrows():
        gene = data[f'{root_gene}_connections'] 
        # some conditionals depending on regulon value ['None', OR list of regulons >= 1] 
        regulon = data['Reg'] 
        if regulon == 'None':
            pass
        else:
            G.add_node(gene, color = 'yellow') 
            for regulon in regulon:
                G.add_node(regulon, color = 'red') 
                G.add_edge(regulon, gene)
            G.add_edge(root_gene,gene)
    colors = [G.nodes[node]['color'] for node in G.nodes]
    fig, ax = plt.subplots(1,1,figsize = (10,10))
    ax.set_title(f'{root_gene} subnetwork.', fontsize = 12, fontweight = 2)
    nx.draw(G,with_labels = True, node_color = [G.nodes[node]['color'] for node in G.nodes], ax = ax)
    nx.write_gexf(G, f'{path}{root_gene}_target_network.gexf')

def remove_sex_genes(adata,species):
    annot = sc.queries.biomart_annotations(
        species,
        ["ensembl_gene_id","external_gene_name","start_position",
         "end_position","chromosome_name"],).set_index("external_gene_name")
    y_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
    x_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "X"])
    
    sex_genes = y_genes.union(x_genes)
    #print(len(sex_genes))
    all_genes = adata.var.index.tolist()
    #print(len(all_genes))
    
    keep_genes = [x for x in all_genes if x not in sex_genes]
    #print(len(keep_genes))
    
    adata_filtered = adata[:,keep_genes] 
    print(f'Removed {len(sex_genes)} sex genes. {len(keep_genes)} remain.')
    return adata_filtered
    
# For specific regulon genes. 
def GO_network_regulon(regulon_object,adata,regulon,direction:str,save:str,attributes:list,path:str,
                       groupby:str='genotype',species:str="mmusculus",selection:str='Wt',
                       age:str='P3',manual_list:list=None,deg:bool=False, p_cutoff:int=0.05,lfc=0.5):
    '''
    t-test overestimated DEG on regulon-target subset genes. Or if no regulon is set, do broad GSEA. 
    Default groupby does DEG analysis on basis of genotype.

    User can input gene list in place of regulon targets by passing a list of genes to 'manual_list' arg. 
    ''' 
    warnings.simplefilter(action='ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', message='Initializing view as actual')
    warnings.filterwarnings('ignore', message='Trying to set attribute \".obs\" of view, copying.')

    ro = regulon_object
    a = adata.copy()
    dlfc = lfc*-1 # inverse for negative log2fc threshold.
    # For user-defined gene list.
    if manual_list !=None:
        genes = manual_list
        a = adata[:, genes]
        if len(age) == 1:
            a = a[a.obs['age']==age[0]] 
        else:
            a = a[a.obs['age'].isin(age)]
        # remove sex genes, and perform gseapy+enrichR GSEA analysis.
        a = remove_sex_genes(a, species) 
        sc.tl.rank_genes_groups(a, groupby,method = 't-test_overestim_var',key_added="t-test_ov")
        gene_set_names = gseapy.get_library_name(organism='Mouse')
        if direction == 'up':
            print(f'Building network from upregulated genes in {selection}')
            glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_min=lfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
            print(f'Glist has {len(glist)} genes')
        if direction == 'down':
            print(f'Building network from downregulated genes in {selection}')
            glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_max=dlfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
            print(f'Glist has {len(glist)} genes')
        enr_res = gseapy.enrichr(gene_list=glist,organism="Mouse",
                                 gene_sets="GO_Biological_Process_2021",
                                 cutoff = 0.01)
        # GSEA conversion into network table.
        nodes, edges = gseapy.enrichment_map(enr_res.res2d,top_term = 300,column='Adjusted P-value')
        network_table = edges[['src_name','targ_name','jaccard_coef','overlap_coef','overlap_genes']]
        network_table.columns = ['source','target','jaccard','overlap','intersection']
        G = nx.from_pandas_edgelist(network_table, source = 'source',
                                target='target',
                                edge_attr = ['jaccard','overlap','intersection']) 
        node_attributes = {node: pval for node, pval in zip(nodes['Term'],nodes['Adjusted P-value'])}

        # P-value of enrichment score will be a node attribute (size,color, etc.) 
        nx.set_node_attributes(G,node_attributes,"pval")
        nx.write_gexf(G, f'{path}{save}.gexf')
        
    else:
         # Create gene subset based on regulon targets
        if regulon !=None: 
            # For looking at a singular regulon. 
            if type(regulon) == str:
                r = regulon+'(+)' # pySCENIC suffix.
                i = find_regulon(ro,r, verbose = False) 
                t = list(ro[i].gene2weight.keys()) 
                a = adata[:, t] # subset by gene list 
                print(a)
                if len(age) == 1:
                    a = a[a.obs['age'] == age[0]]
                else:
                    a = a[a.obs['age'].isin(age)]
                # Remove sex-chromosome genes via biomart query. 
                a = remove_sex_genes(a,species)
                # Rank genes given group_by condition using scanpy.
                sc.tl.rank_genes_groups(a, groupby,method = 't-test_overestim_var',key_added="t-test_ov")
                gene_set_names = gseapy.get_library_name(organism='Mouse')
                if direction == 'up':
                    print(f'Building network from upregulated genes in {selection}')
                    glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_min=lfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                    print(f'Glist has {len(glist)} genes')
                if direction == 'down':
                    print(f'Building network from downregulated genes in {selection}')
                    glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_max=dlfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                    print(f'Glist has {len(glist)} genes')
                enr_res = gseapy.enrichr(gene_list=glist,organism="Mouse",
                                         gene_sets="GO_Biological_Process_2021",
                                         cutoff = 0.01)
                
                nodes, edges = gseapy.enrichment_map(enr_res.res2d,top_term = 300,column='Adjusted P-value')
                network_table = edges[['src_name','targ_name','jaccard_coef','overlap_coef','overlap_genes']]
                network_table.columns = ['source','target','jaccard','overlap','intersection']
                G = nx.from_pandas_edgelist(network_table, source = 'source',
                                        target='target',
                                        edge_attr = ['jaccard','overlap','intersection']) 
                node_attributes = {node: pval for node, pval in zip(nodes['Term'],nodes['Adjusted P-value'])}
                
                nx.set_node_attributes(G,node_attributes,"pval")
                nx.write_gexf(G, f'{path}{save}.gexf')
            if type(regulon) == list: 
                #If you are using multiple regulons' target genes. 
                targets = [] 
                # Take a union of regulon targets in list. 
                for reg in regulon:
                    r = reg+'(+)' # pySCENIC suffix.
                    i = find_regulon(ro,r, verbose = False)
                    t = list(ro[i].gene2weight.keys())
                    targets.extend(t) 
                t = list(set(targets))
                a = adata[:,t]
                if len(age) == 1:
                    a = a[a.obs['age'] == age[0]]
                else:
                    a = a[a.obs['age'].isin(age)]
                # Remove sex genes via biomart query. 
                a = remove_sex_genes(a, species)
                # Rank genes given group_by condition using scanpy. 
                sc.tl.rank_genes_groups(a, groupby,method = 't-test_overestim_var',key_added="t-test_ov")
                gene_set_names = gseapy.get_library_name(organism='Mouse')
                if direction == 'up':
                    print(f'Building network from upregulated genes in {selection} cells')
                    glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_min=lfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                    print(f'Glist has {len(glist)} genes')
                if direction == 'down':
                    print(f'Building network from downregulated genes in {selection} cells')
                    glist = sc.get.rank_genes_groups_df(a, group = selection,key="t-test_ov",log2fc_max=dlfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                    print(f'Glist has {len(glist)} genes')
                enr_res = gseapy.enrichr(gene_list=glist,organism="Mouse",
                                         gene_sets="GO_Biological_Process_2021",
                                         cutoff = 0.01)
                
                nodes, edges = gseapy.enrichment_map(enr_res.res2d,top_term = 300,column='Adjusted P-value')
                network_table = edges[['src_name','targ_name','jaccard_coef','overlap_coef','overlap_genes']]
                network_table.columns = ['source','target','jaccard','overlap','intersection']
                G = nx.from_pandas_edgelist(network_table, source = 'source',
                                        target='target',
                                        edge_attr = ['jaccard','overlap','intersection']) 
                
                node_attributes = {node: pval for node, pval in zip(nodes['Term'],nodes['Adjusted P-value'])}
                nx.set_node_attributes(G,node_attributes,"pval")
                nx.write_gexf(G, f'{path}{save}.gexf')
        else:
            # If you do not have a gene list or regulon list. perform DEG analysis. 
            print('Performing DEG analysis prior to GSEA')
            if len(age) == 1:
                a = a[a.obs['age'] == age[0]]
            else:
                a = a[a.obs['age'].isin(age)]
            a = remove_sex_genes(a,species)
            sc.tl.rank_genes_groups(a, groupby, method ='t-test_overestim_var',key_added = "t-test_ov") 
            gene_set_names = gseapy.get_library_name(organism='Mouse')
            if direction == 'up':
                print(f'Building network from upregulated genes in {selection} cells')
                glist = sc.get.rank_genes_groups_df(a,group = selection, key = "t-test_ov",log2fc_min=lfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                print(f'Glist has {len(glist)} genes')
            if direction == 'down':
                print(f'Building network from downregulated genes in {selection} cells')
                glist = sc.get.rank_genes_groups_df(a,group = selection, key = "t-test_ov",log2fc_max=dlfc,pval_cutoff=p_cutoff)['names'].squeeze().str.strip().tolist()
                print(f'Glist has {len(glist)} genes')
            enr_res = gseapy.enrichr(gene_list=glist,organism='Mouse',
                                     gene_sets = 'GO_Biological_Process_2021',
                                     cutoff = 0.01)
            nodes,edges = gseapy.enrichment_map(enr_res.res2d,top_term=300,column='Adjusted P-value')
            network_table = edges[['src_name','targ_name','jaccard_coef','overlap_coef','overlap_genes']]
            network_table.columns = ['source','target','jaccard','overlap','intersection']
            G = nx.from_pandas_edgelist(network_table,source = 'source',
                                        target = 'target',
                                        edge_attr = ['jaccard','overlap','intersection'])
            node_attributes = {node: pval for node, pval in zip(nodes['Term'],nodes['Adjusted P-value'])}
            nx.set_node_attributes(G,node_attributes,"pval")
    
    # Add node attribute for intersection to regulon as proportion of targets 
    overlap_attributes = {}
    for collection in attributes:
        overlap_attributes[collection] = {} 
        print(f'Calculating intersection for {collection} targets.')
        i = find_regulon(regulon_object,collection,verbose=False)
        targets = list(regulon_object[i].gene2weight.keys())
        targets = map(lambda x: x.upper(), targets)
        targets = set(targets)
        for row, data in enr_res.res2d.iterrows():
            term = data['Term']
            overlap = set(data['Genes'].split(';'))
            intersection = overlap.intersection(targets)
            repRatio = len(intersection) / len(targets)
            overlap_attributes[collection][term] = repRatio
        node_attributes = {node: ratio for node, ratio in zip(overlap_attributes[collection].keys(), overlap_attributes[collection].values())}
        nx.set_node_attributes(G,node_attributes,f'{collection}_ratio')
    nx.write_gexf(G, f'{path}{save}.gexf')
            
            
        
    n_nodes = len(nodes)
    n_edges = len(network_table)
    print(f'Network of {n_nodes} nodes and {n_edges} edges generated at path.') 
    return enr_res