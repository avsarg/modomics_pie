# -*- coding: utf-8 -*-
"""
@author: avsarg
"""
import pandas as pd
import matplotlib.pyplot as plt

def modomics_pie():
    RNA = pd.read_csv('Modomics - A Database of RNA Modifications.csv', sep=',')
    
    data = RNA.copy()
    
    #redesign the modification names
    for i in range(len(data)):
        if data['Reaction'][i] =='A:m6A':data['Reaction'][i] = 'A<->m6A'
        if data['Reaction'][i] =='xX:Xm':data['Reaction'][i] = 'N->Nm'
        else: data['Reaction'][i]=data['Reaction'][i].replace(':', '->')
    
    
    # info for the number of occurances of each diseases
    from collections import Counter
    occurrences = Counter(data['Disease Name'])
    occurrences = pd.DataFrame.from_dict(occurrences, orient='index')
    occurrences.columns = ['Occurr']
    
    # Define cancer group diseases
    cancs = [j for j in occurrences.index.to_list() if ('ancer' in j) == True]
    Leuks = [j for j in occurrences.index.to_list() if ('Leuk' in j) == True]
    omas = [j for j in occurrences.index.to_list() if ('oma' in j) == True]
    
    all_cancs = cancs + Leuks + omas
    
    # group the diseases as cancerous or non-cancerous
    data['Grp'] = ''
    for i in data['Disease Name']:
        if i in all_cancs:
            data.loc[data['Disease Name']==i,'Grp'] = 'cancerous'
        else:
            data.loc[data['Disease Name']==i,'Grp'] = 'non-cancerous'
    
    # a list for the diseasees first cancerous (28) then others (40)
    for i in occurrences.index.to_list():
        if (i in all_cancs) == False:
            all_cancs.append(i)
    
    # create a df the the order of the diseases are the same with all_cancs list
    occurr_allDis = pd.DataFrame(columns=['sorterd_dis', 'occurrence'])
    for i in range(len(all_cancs)):
        occurr_allDis.loc[i,'sorterd_dis'] = all_cancs[i]
        occurr_allDis.loc[i,'occurrence'] = occurrences.loc[all_cancs[i]].values[0]
    
    #make the disease names index
    occurr_allDis = occurr_allDis.set_index('sorterd_dis')
    
    
    # For the occurances of the reactions
    occurr_rxn = Counter(data['Reaction'])
    occurr_rxn = pd.DataFrame.from_dict(occurr_rxn, orient='index')
    occurr_rxn.columns = ['Occurr_rxn']
    
    
    # create occurrences table for each rxn and the diseases (7,68)
    df = data[['Reaction', 'Grp']]
    # Group boroughs and crimes and count numbers
    df_grouped = df.groupby(['Reaction', 'Grp']).size().reset_index()
    
    # Create a pivot table
    table = pd.pivot_table(df_grouped, index=['Reaction'], columns=['Grp'])
    # Delete columns with Nan values and convert to integers
    table = table.fillna(0)
    table
    
    table.columns = [table.columns.to_flat_index()[i][1] for i in range(len(table.columns)) ]
    table
    
    # rearrange the order of the index
    table = table.reindex(index = occurr_rxn.index)
    table
    
    # Define the groups: I have 7 maingroups (at the inner ring) and 68 subgroups (at the outer ring)
    group_names= occurr_rxn.index.to_list()
    group_size= occurr_rxn.to_numpy().flatten()
    
    subgroup_names = [i+" : "+j for i in group_names for j in table.columns.to_list()]
    #subgroup_size=[84,17,3,7,16,12,1,7,4,6,5,2,8,3]
    subgroup_size = []
    for i in range(len(table)):
        vals = table.iloc[i,:].to_list()
        #size.append(vals)
        subgroup_size = subgroup_size + vals
    
    # prepare the colors for outer
    colr = ['lightsteelblue', 'plum']
    sub_color = colr*7
    
    # prepare the colors for inner
    colr = ['#cc0066', '#4d0099','#40507d', '#006699','#00994d', '#6b6b47','#993300']
    
    # First Ring (outside)
    fig, ax = plt.subplots(figsize=(18,15))
    ax.axis('equal')
    mypie, _ = ax.pie(subgroup_size, radius=1.3, labels=subgroup_names,
                      colors=sub_color, labeldistance=None,
                      textprops={'fontsize': 10, 'weight':'bold'},
                      wedgeprops={'edgecolor': 'white', 'linewidth': 2})
    plt.setp( mypie, width=0.5, edgecolor='white')
    plt.margins(0,0)
    
    # Second Ring (Inside)
    mypie2, _ = ax.pie(group_size, radius=1.3-0.3,
                       labels=group_names, labeldistance=0.53, colors=colr, rotatelabels=True,
                       textprops={'fontsize': 16, 'weight':'bold', 'color': 'white'},
                       wedgeprops={'edgecolor': 'white', 'linewidth': 2})
    plt.setp( mypie2, width=0.5, edgecolor='white')
    plt.margins(0,0)
    
    '''# for labels
    handles, labels = ax.get_legend_handles_labels()
    subgroup_names_legs=['Cancerous','Non-cancerous']
    ax.legend(handles[:2], subgroup_names_legs, loc=(0.95, 0.6), fontsize=20)'''
    
    plt.tight_layout()
    fig.savefig('Modomics_NewDiseases_piechart.pdf', format='pdf', dpi=1200,bbox_inches='tight', transparent=True)
    plt.show()
    
    return(fig)


