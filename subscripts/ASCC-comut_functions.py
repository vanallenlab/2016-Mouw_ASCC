# Brendan Reardon
# Van Allen Laboratory
# Dana-Farber Cancer Institute, Harvard Medical School
# The Broad Institute of MIT & Harvard
# 10 April 2016

# Genomic Evolution After Chemoradiotherapy in Anal Squamous Cell Carcinoma
# CoMut Functions

## ============================================================================================
# A function to sort sample names for each individual by their tumor type
# this will be utilized to stack our CoMut plots in Figure 1

def populate_samples_by_type(info_dataframe, clinical_outcome):
    outcome_dataframe = info_dataframe[info_dataframe.ix[:,'Clinical Outcome'] == clinical_outcome]
    info_dataframe_individuals = outcome_dataframe.ix[:,'Individual ID']
    ids_dataframe = pd.DataFrame([], columns = info_dataframe_individuals.unique().tolist(), index = ['Primary', 'Recurrent'])
    
    for individual in ids_dataframe.columns.tolist():
        tmp = info_dataframe[info_dataframe.ix[:,'Individual ID'] == individual]
        for i in tmp.index.tolist():
            tmp_type = tmp.ix[i, 'Tumor Type']
            ids_dataframe.ix[tmp_type, individual] = tmp.ix[i, 'Tumor Sample ID']

    return(ids_dataframe)

## ============================================================================================
# A function to populate the dataframes above with numbers representing each variant type

def populate(list_samples, df_genes, cohort_maf):
    # We create a dataframe
    dataframe = pd.DataFrame(index = range(0, len(df_genes)), columns = [['gene'] + list_samples]).fillna(0)
    dataframe.ix[:, 'gene'] = df_genes['gene']

    # We subset for only our samples of interest and respective mutations
    variants = cohort_maf[cohort_maf['Tumor_Sample_Barcode'].isin(list_samples)]
    
    # And populate the appropriate element 
    for i in variants.index.tolist():
        r = dataframe[dataframe['gene'] == variants.ix[i, 'Hugo_Symbol']].index.tolist()[0]
        c = variants.ix[i, 'Tumor_Sample_Barcode']
        
        if variants.ix[i, 'Variant_Classification'] == 'Silent':
            dataframe.loc[r, c] = 1
        elif variants.ix[i, 'Variant_Classification'] == 'Missense_Mutation':
            dataframe.loc[r, c] = 2
        elif variants.ix[i, 'Variant_Classification'] == 'Splice_Site':
            dataframe.loc[r, c] = 3        
        elif variants.ix[i, 'Variant_Classification'] == 'Nonsense_Mutation':
            dataframe.loc[r, c] = 4
        elif variants.ix[i, 'Variant_Classification'] == 'Frame_Shift_Del':
            dataframe.loc[r, c] = 5        
        elif variants.ix[i, 'Variant_Classification'] == 'Frame_Shift_Ins':
            dataframe.loc[r, c] = 5
        elif variants.ix[i, 'Variant_Classification'] == 'In_Frame_Del':
            dataframe.loc[r, c] = 6        
        elif variants.ix[i, 'Variant_Classification'] == "3'UTR":
            dataframe.loc[r, c] = 7
        elif variants.ix[i, 'Variant_Classification'] == 'Intron':
            dataframe.loc[r, c] = 8

    return(dataframe)

## ============================================================================================
# A function used to reorder dataframes horizontally and vertically based on the gene mutation 
# frequency across samples
def format_freq_sort(dataframe, frequency, number_of_genes):
    # Sort dataframe across genes depending on frequency across all cohorts
    dataframe_sorted = pd.DataFrame(dataframe, index = frequency.index)
    # Reset index
    dataframe_sorted.index = range(0, len(dataframe_sorted))
    # Sort across samples for aesthic purposes
    dataframe_sorted = dataframe_sorted.T.sort_values(by = range(number_of_genes), ascending = [0] * number_of_genes).T
    dataframe_sorted.index = range(0, len(dataframe_sorted))

    return(dataframe_sorted)

## ============================================================================================
# A function to populate a dataframe with the recurrence location
# Used for figure 1, clinical annotation 1: primary samples - recurrent status

def format_resp(info_cohort, input_dataframe):
    output_dataframe = pd.DataFrame([], columns = input_dataframe.columns.tolist(), index = ['Recurrence'])

    for sample in output_dataframe.columns.tolist():
        if (isinstance(sample, str)):
            index = info_cohort[info_cohort.ix[:, 'Tumor Sample ID'] == sample].index.tolist()[0]
            status = info_cohort.ix[index, 'Clinical Outcome']
            if status == 'Response':
                num_loc = 1
            else:
                num_loc = 0
            output_dataframe.ix['Recurrence', sample] = num_loc
        else:
            output_dataframe.ix['Recurrence', sample] = 10
            
    return(output_dataframe)

## ============================================================================================
# A function to populate a dataframe with the recurrence location
# Used for figure 1, clinical annotation 1: recurrent samples - recurrent location

def format_loc(info_cohort, input_dataframe):
    output_dataframe = pd.DataFrame([], columns = input_dataframe.columns.tolist(), index = ['Recurrence Location'])

    for sample in output_dataframe.columns.tolist():
        if (isinstance(sample, str)):
            index = info_cohort[info_cohort.ix[:, 'Tumor Sample ID'] == sample].index.tolist()[0]
            loc = info_cohort.ix[index, 'Recurrence Distance']
            if loc == 'distant':
                num_loc = 1
            else:
                num_loc = 0
            output_dataframe.ix['Recurrence Location', sample] = num_loc
        else:
            output_dataframe.ix['Recurrence Location', sample] = 10

    return(output_dataframe)

## ============================================================================================
# A function used to populate a dataframe with the hpv status
# Used for supplementary figure 4, clinical annotation 1: primary samples - cohort status

def format_cohort(input_dataframe_0, input_dataframe_1):
    zero = np.zeros(len(input_dataframe_0.columns.tolist()))
    one = np.ones(len(input_dataframe_1.columns.tolist()))
    np_output = np.concatenate((zero, one), axis = 0)

    output_dataframe = pd.DataFrame(np_output).T
    output_dataframe.columns = input_dataframe_0.columns.tolist() + input_dataframe_1.columns.tolist()
    output_dataframe = output_dataframe = output_dataframe.rename(index = {0: 'Cohort Status'})

    return(output_dataframe)

## ============================================================================================
# A function used to populate a dataframe with the hpv status
# Used for figure 1, clinical annotation 2: all samples - hpv status

def format_hpv(info_cohort, input_dataframe):
    output_dataframe = pd.DataFrame([], columns = input_dataframe.columns.tolist(), index = ['HPV Status'])

    for sample in output_dataframe.columns.tolist():
        if (isinstance(sample, str)):
            index = info_cohort[info_cohort.ix[:, 'Tumor Sample ID'] == sample].index.tolist()[0]
            status = info_cohort.ix[index, 'HPV Status']
            if status == 'Negative':
                num_status = 0
            elif status == 'Positive':
                num_status = 1
            else: 
                num_status = 2 # For the case of unknown HPV status (n = 1)
        else:
            num_status = 10

        output_dataframe.ix[0, sample] = num_status

    return(output_dataframe)

## ============================================================================================
# A function used to populate a dataframe with the 3q status from GISTIC
# Used for figure 1, genomic annotation 1: all samples - 3q CNV

def format_gistic_3q(cohort_3q, input_dataframe):
    output_dataframe = pd.DataFrame([], columns = input_dataframe.columns.tolist(), index = ['3q Status'])

    for sample in output_dataframe.columns.tolist():
        if (isinstance(sample, str)):
            value = cohort_3q[sample]
            if value > 0:   # Gain
                output_dataframe.ix['3q Status', sample] = 1 
            elif value < 0: # Loss
                output_dataframe.ix['3q Status', sample] = -1
            else:
                output_dataframe.ix['3q Status', sample] = 0
        else:
            output_dataframe.ix['3q Status', sample] = 10

    return(output_dataframe)

## ============================================================================================
# A function used calculate the mutational burden of the pilot cohort
# Used in figure 1

def format_pilot_load(mutsig_rates, input_dataframe):
    rows = ['Total', 'Nonsyn', 'Syn']
    output_dataframe = pd.DataFrame([], columns = input_dataframe.columns.tolist(), index = rows)

    for sample in output_dataframe.columns.tolist():
        if (isinstance(sample, str)):
            index = mutsig_rates[mutsig_rates.ix[:,'name'] == sample].index.tolist()[0]
            # We extract our values and convert from Mutations/base to Mutations/Megabase
            output_dataframe.ix['Total', sample] = mutsig_rates.ix[index, 'rate_tot'] *10*10*10*10*10*10
            output_dataframe.ix['Nonsyn', sample] = mutsig_rates.ix[index, 'rate_non'] *10*10*10*10*10*10
            output_dataframe.ix['Syn', :] = output_dataframe.ix['Total', :] - output_dataframe.ix['Nonsyn', :]
        else:
            output_dataframe.ix['Total', sample] = 0
            output_dataframe.ix['Nonsyn', sample] = 0
            output_dataframe.ix['Syn', :] = 0

    return(output_dataframe)

## ============================================================================================
# A function used calculate an estimated non-synonymous mutational burden from unmatched tumor sample MAFs
# Used in figure 1, supplementary figure 4,  for the extension cohort

def format_extension_load(info_cohort, mutations_cohort, input_dataframe):
    samples = input_dataframe.columns.tolist()

    # Create an output dataframe
    tmp_dataframe = pd.DataFrame([], columns = ['Tumor Sample ID', 'Exome Somatic Bases Covered', '# Mutations', '# Mutations / Mb'],
        index = range(0, len(samples)))
    tmp_dataframe['Tumor Sample ID'] = samples

    for sample in samples:
        if (isinstance(sample, str)):
            # Define our indicies
            index = info_cohort[info_cohort['Tumor Sample ID'] == sample].index.tolist()[0]
            index_tmp = tmp_dataframe[tmp_dataframe['Tumor Sample ID'] == sample].index.tolist()[0]

            # Exome Somatic Bases
            tmp_dataframe.ix[index_tmp, 'Exome Somatic Bases Covered'] = info_cohort.ix[index, 'Exome Somatic Bases Covered']

            # # Mutations
            df_tmp = mutations_cohort[mutations_cohort['Tumor_Sample_Barcode'] == sample]
            tmp_dataframe.ix[index_tmp, '# Mutations'] = len(df_tmp)
            tmp_dataframe.ix[index_tmp, '# Mutations / Mb'] = tmp_dataframe.ix[index_tmp, '# Mutations'] / tmp_dataframe.ix[index_tmp, 'Exome Somatic Bases Covered']
            tmp_dataframe.ix[index_tmp, '# Mutations / Mb'] = tmp_dataframe.ix[index_tmp, '# Mutations / Mb']*10*10*10*10*10*10
        else:
            # Define our indicies
            index_tmp = samples.index(sample)
            tmp_dataframe.ix[index_tmp, '# Mutations / Mb'] = 0

    output_dataframe = tmp_dataframe.ix[:,'# Mutations / Mb'].to_frame('Nonsyn')
    output_dataframe.index = samples
    output_dataframe = output_dataframe.T

    output_dataframe.ix['Syn', :] = 0

    return(output_dataframe)

## ============================================================================================
# A function used to create a CoMut plot with matplotlib
# Used in figure 1, supplementary figure 4

def plot_comut(gridspec_y, gridspec_x, dataframe, display_label, toggle_label, label_style, color_cmap, color_norm, line_color, 
    sharex = None, sharey = None):
    # Create gridspec location for comut
    ax = plt.subplot(gs[gridspec_y, gridspec_x], sharex = sharex, sharey = sharey)

    # Grab the number of samples and length of dataframe
    nsamples = len(dataframe.columns.tolist())
    length_dataframe = len(dataframe)


    # Create lines across the grid
    for y in range(-1, length_dataframe, 1):
            plt.plot(range(-1, nsamples + 1, 1), [y + 0.5] * len(range(-1, nsamples + 1, 1)), '-', lw = 2, color = line_color, alpha = 1)
    for x in range(-1, nsamples, 1):
            plt.plot([x + 0.5] * len(range(-1, length_dataframe + 1)), range(-1, length_dataframe + 1), '-', lw = 2, color = line_color, alpha = 1)
    # Remove the plot frame lines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)


    # Specify toggle for if gene list should be on or off
    if toggle_label == 1:
        label_left = 'on'
    else: 
        label_left = 'off'

    # Remove the tick marks
    plt.tick_params(axis = 'both', which = 'both',
        top = 'off', bottom = 'off', right = 'off', left = 'off',
        labeltop = 'off', labelbottom = 'off', labelright = 'off', labelleft = label_left) # Will toggle to show gene

    # Resize tick marks
    plt.xticks(range(0, nsamples, 1))
    plt.yticks(range(0, length_dataframe, 1), fontsize = 13)

    # Replace y ticks with string values (gene name)
    ax.set_yticklabels(display_label, style = label_style)

    cax = ax.imshow(dataframe.astype(int), interpolation = 'none', aspect = 'auto', cmap  = color_cmap, norm = color_norm)

    return(ax, cax)

def plot_burden(gridspec_y, gridspec_x, dataframe, toggle_label, colors, ylow, ymax,
    display_label = None, sharex = None, sharey = None, label_bottom = 'off'):
    # Create gridspec location for bar graph
    ax = plt.subplot(gs[gridspec_y, gridspec_x], sharex = sharex, sharey = sharey)

    # Remove the plot frame lines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Specify toggle for if gene list should be on or off
    if toggle_label == 1:
        label_left = 'on'
    else: 
        label_left = 'off'

    # Remove the tick marks
    plt.tick_params(axis = 'both', which = 'both', 
                top = 'off', bottom = 'off', right = 'off', left = 'off',
                labeltop = 'off', labelbottom = label_bottom, labelright = 'off', labelleft = label_left)

    # Nonsyn burden
    axbar = ax.bar(range(0, len(dataframe.columns.tolist())), dataframe.ix['Nonsyn', :],
        width = 0.85, align = 'center')
    # Syn burden
    axbar2 = ax.bar(range(0, len(dataframe.columns.tolist())), dataframe.ix['Syn', :],
        width = 0.85, align = 'center', bottom = dataframe.ix['Nonsyn', :])

    # Color lines across the grid
    for y in range(ylow, ymax, 10):    
        plt.plot(range(-1, len(dataframe.columns.tolist()) + 1), [y] * (len(range(-1, len(dataframe.columns.tolist()))) + 1), 
             "--", lw=0.5, color="black", alpha=0.6)

    plt.xlim(-0.5, len(dataframe.columns.tolist()) - 0.5)
    plt.ylim(0, ymax)

    for i in range(0, len(dataframe.columns.tolist())):
        axbar[i].set_color(colors_bar[0])
        axbar[i].set_edgecolor('black')
        axbar2[i].set_color(colors_bar[1])
        axbar2[i].set_edgecolor('black')

    # Set ylabel
    if display_label != None:
        ax.set_ylabel(display_label, fontsize = 12, rotation = 0, labelpad = 50)

    return(ax)


def create_legend(gridspec_y, gridspec_x, labels, colors, bbox_x, bbox_y, ncols): 
    # Create gridspec location for the legend
    ax = plt.subplot(gs[gridspec_y, gridspec_x])

    # Remove the plot frame lines
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tick_params(axis = 'both', which = 'both', 
                top = 'off', bottom = 'off', right = 'off', left = 'off',
                labeltop = 'off', labelbottom = 'off', labelright = 'off', labelleft = 'off')

    for label_ in labels:
        count = labels.index(label_)
        ax.bar(count, 0, label = label_, color = colors[count])

    plt.ylim(0)
    legend = ax.legend(bbox_to_anchor = (bbox_x, bbox_y), ncol = ncols, handlelength = 3, fontsize = 14, 
        columnspacing = 0.5, frameon = False)

    legend.get_title().set_fontsize('14')

def fillna_columns(input_dataframe):
    columns_input = input_dataframe.columns.tolist()
    
    # Break up columns
    columns_samples = [column_ for column_ in columns_input if column_ == column_]
    columns_nan = [column_ for column_ in columns_input if column_ != column_]

    # Create dataframe of 10s for nan column titles
    df_nan = pd.DataFrame(index = range(0, len(input_dataframe)), columns = range(0, len(columns_nan)))
    df_nan = df_nan.fillna(10)
    
    # Combine with previous dataframe
    output_dataframe = pd.concat([input_dataframe.ix[:,columns_samples], df_nan], axis = 1)
    
    return(output_dataframe)

def list_remove_nan(input_list):
    list_ = [x for x in input_list if x == x]

    return(list_)


