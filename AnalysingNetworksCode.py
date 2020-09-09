# -*- coding: utf-8 -*-
"""
@authors: Shaurya Dev Singh
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import collections
import random
import EoN as eon
from scipy.stats import kendalltau
from math import ceil, floor
from functools import partial
from collections import namedtuple
from functools import wraps

# Generating synthetic versions of Colorado Springs Network with Configuration Model graphs
edgelist = pd.read_csv("edgelist.truecolsprings.csv")
edgelist.head()
edgelist.shape

# Initialising empty graph object (undirected)
G = nx.Graph()


# loop to create edgelist in proper format for networkx
edgelist_format = []
for index, row in edgelist.iterrows(): # to get row indices of the data frame
    inside = []
    inside.append(edgelist.at[index,'V1']) 
    inside.append(edgelist.at[index,'V2'])
    edgelist_format.append(inside) # appending pairs to main list iteratively
edgelist_format[:10]
    # inside only gets stored for the most recent iteration 

len(edgelist_format)
# length is corect and last entry is correct - verified with data frame/csv file

# Adding edges to the graph
G.add_edges_from(edgelist_format)
nx.info(G) # getting the correct number of edges in the graph object (number of nodes also coresponds to the max value
           # from the csv file)

# degree sequence sorted in descending order
degree_sequence = sorted([d for n, d in G.degree()], reverse=True)   
degree_sequence[:50]

# Sampling 1000 nodes from the degree sequence
degree_sequence_sample = random.sample(degree_sequence, 1000)
degree_sequence_sample.sort(reverse = True)

# degree distribution for the true network
degreeCount1 = collections.Counter(degree_sequence)           
degr, cont = zip(*degreeCount1.items())
plt.bar(degr, cont, width=0.80, color='b')
plt.title("Degree Histogram")
plt.ylabel("Count")
plt.xlabel("Degree")
plt.yscale('log')
plt.show()
# run these lines together to get desired graph


# generating multiple networks initially using undirected configuration models
graphs = []
while len(graphs) < 500: # number of copies required
    graphs.append(nx.configuration_model(degree_sequence_sample))    
print(nx.info(graphs[0]), nx.info(graphs[3])) # graph structure will be same but there will be randomness in nodes and edges

# removing parallel edges
nx.info(graphs[0])
graphs2 = list(map(nx.Graph, graphs))
print(nx.info(graphs[4]), nx.info(graphs2[4])) # second one should have a lower number edges

graphs3 = []
for g in graphs2:
    g_new = g.copy()
    g_new.remove_edges_from(nx.selfloop_edges(g_new))
    graphs3.append(g_new)
    

print(nx.info(graphs[4]), nx.info(graphs2[4]), nx.info(graphs3[4]))
# NOTE: multiple edges and self loops are allowed in 'multigraph' type objects, not in 'graph' type objects

# Getting a histogram for the degree distribution of one configuration model random graph
degree_sequence_conf = sorted([d for n, d in graphs3[1].degree()], reverse=True)  # degree sequence of the configuration model
degreeCount = collections.Counter(degree_sequence_conf)                           # random graph
deg, cnt = zip(*degreeCount.items())
plt.bar(deg, cnt, width=0.80, color='b')
plt.title("Degree Histogram")
plt.ylabel("Count")
plt.xlabel("Degree")
plt.yscale('log')
plt.show()
# identical as true network

# removing nodes at random - 10% to 90% with breaks of 10

# obtaining the node list - will not change based on the graphs, so we take it from any one
node_list = list(graphs3[0].nodes)
len(node_list)
# nx.info(G1)

# creating a loop to generate the sample nodes for us
perc = [.10,.20,.30,.40,.50,.60,.70,.80,.90]
random_sample_list = [] 
for p in perc:
    sample = round(p*len(graphs3[0].nodes))
    interior_list = []
    while len(interior_list) < 500: # number of copies required
        random_sample = random.sample(node_list, sample)
        interior_list.append(random_sample)
    random_sample_list.append(interior_list)


# some checks - to see if code worked properly
len(random_sample_list)    # 9 sub lists created for different missing percentages
len(random_sample_list[0]) # each missing level has 10 random samples for the 10 graphs
# the second one shold be more - 10% and 20% cases
print(len(random_sample_list[0][1]), len(random_sample_list[1][4])) 


# getting the graph object list after removals
graph_list = []
for i in range(len(random_sample_list)): # keeps percentage removal counter
    interior_graphs = []    
    for g,j in zip(graphs3, range(len(graphs3))): # zip keeps graph object and graph number counters parallel
        g_new = g.copy() # object explicitly extracted - required to make it work as it is stored explicitly within the loop
        # does not work on only g as it is not updated in the list, with the line below
        g_new.remove_nodes_from(random_sample_list[i][j])
        interior_graphs.append(g_new)
    graph_list.append(interior_graphs) # after inner for loop is finished, append filled interior vector outside


print(nx.info(graph_list[0][3]), nx.info(graph_list[1][4])) # gives different number of final nodes
print(len(graph_list), len(graph_list[3])) # 9 percentages and 10 each


# adding base graphs as first entries in list
len(graph_list)
graph_list2 = [graphs3] + graph_list
len(graph_list2)
print(nx.info(graph_list2[0][3]), nx.info(graph_list2[1][4])) # no removals, 10% removals, so on...

# viewing info of all base graphs and 900 removals
result = list(map(nx.info, graph_list2[0]))
result2 = list(map(nx.info, graph_list2[9]))
print(result, result2)

########################################################################################################

# IMMUNISATION ROBUSTNESS

# Getting node removal/immunsation lists for all graphs with their respective percentages of removed nodes
q = 0.10 # setting the parameter for percentage immunised


# Random Immunisation - basis for comparison with each strategy


immunised_graph_list_random = []
for j in range(len(graph_list2)): # for all percentage brackets
    inner_immunised_list = []
    for gra in graph_list2[j]: # for all graphs in one percentage bracket
        g1 = gra.copy()
        node_list2 = list(g1.nodes)
        sampled_nodes = random.sample(node_list2, round(len(node_list2)*q))
        g1.remove_nodes_from(sampled_nodes)
        inner_immunised_list.append(g1)
    immunised_graph_list_random.append(inner_immunised_list)


# code checks
print(nx.info(immunised_graph_list_random[0][2]), nx.info(immunised_graph_list_random[0][4]))
# graphs at 0% missing nodes [0] will have 1000-100 = 900 nodes as 10% are removed at random
print(nx.info(immunised_graph_list_random[1][3]), nx.info(immunised_graph_list_random[1][2]))
print(nx.info(immunised_graph_list_random[1][3]), nx.info(immunised_graph_list_random[2][2]))
# graphs at 10% missing level [1] have 1000-100 = 900 nodes, then 10% are removed at random => 
# 900 - 90 = 810
# graphs at 20% missing level [1] have 1000-200 = 800 nodes, then 10% are removed at random => 
# 800 - 80 = 720, then so on.......


# writing a for loop to get simulations (mean and sd stored) for each graph in our list 

nodes = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100] # for plot-divide by respective population

# 100 simulations for each graph
outbreak_list_random = []
for j, n in zip(range(len(immunised_graph_list_random)), nodes): # for all percentage brackets
    mean_outbreak_list = []
    for g in immunised_graph_list_random[j]: # for all graphs in one percentage bracket
        inner_outbreaks = []
        while len(inner_outbreaks) < 2000:  # 2000 simulations for each graph
            sim = eon.Gillespie_SIR(g, 0.95, 1)
            inner_outbreaks.append((sim[1][0] - sim[1][-1])/n) # S(T) - S(0) (Outbreak Size definition)
        mean_outbreak_list.append(inner_outbreaks)
    outbreak_list_random.append(mean_outbreak_list)


# (1) FOR DEGREE CENTRALITY - just change the first line in the inner loop for others

final_node_immunise_list = []
for j in range(len(graph_list2)): # for all percentage brackets
    inner_immunise_list = []
    for i in range(len(graph_list2[j])): # for all graphs in one percentage bracket
        dc = nx.degree_centrality(graph_list2[j][i])
        dc = pd.DataFrame([dc.keys(), dc.values()]).T
        dc.columns = ['node', 'centrality value'] 
        dc = dc.sort_values(by = 'centrality value', ascending = False)
        dc['node'] = dc['node'].astype(int)
        nodes_to_remove = list(dc['node'].head(round(int(len(dc['node']))*q)))
        inner_immunise_list.append(nodes_to_remove)
    final_node_immunise_list.append(inner_immunise_list) 

# number of immunised nodes (top 10%) should decrease as we go up in missing node graphs because the
# total number of nodes in those grpahs goes down as well


# removing the immunised nodes for the respective graphs
immunised_graph_list = []
for j in range(len(graph_list2)): # for all percentage brackets
    inner_immunised_list = []
    for gra, i in zip(graph_list2[j], range(len(graph_list2[j]))): # for all graphs in one percentage bracket
        g1 = gra.copy()
        g1.remove_nodes_from(final_node_immunise_list[j][i])
        inner_immunised_list.append(g1)
    immunised_graph_list.append(inner_immunised_list)


print(nx.info(immunised_graph_list[0][3]), nx.info(immunised_graph_list[0][4]))
print(nx.info(immunised_graph_list[1][3]), nx.info(immunised_graph_list[2][4]))
# graphs at 10% missing nodes [1] will have 4111-411 = 3700 nodes and then we subtract 
# 10% top nodes that are immunised => 3700 - 370 = 3330 ; etc...


# SIMULATIONS - ROSENBLATT PART - Degree centrality

# writing a for loop to get 100 simulations (mean and sd stored) for each graph in our list 
# take the mean and sd of the 100 first and then take grand mean across 10

nodes = [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100] # for plot-divide by respective population

# 100 simulations for each graph
outbreak_list = []
for j, n in zip(range(len(immunised_graph_list)), nodes): # for all percentage brackets
    mean_outbreak_list = []
    for g in immunised_graph_list[j]: # for all graphs in one percentage bracket
        inner_outbreaks = []
        while len(inner_outbreaks) < 2000:  # 2000 simulations for each graph
            sim = eon.Gillespie_SIR(g, 0.95, 1)
            inner_outbreaks.append((sim[1][0] - sim[1][-1])/n) # S(T) - S(0) (Outbreak Size) (L&F definition)
        mean_outbreak_list.append(inner_outbreaks)
    outbreak_list.append(mean_outbreak_list)


## PLOTTING
outbreak_variation =  np.subtract(outbreak_list_random, outbreak_list)

# putting all the simulations across all graphs into one list, thus removing each graph wise sublist
outbreak_variation_merged = outbreak_variation.reshape(10,1000000) # 10 percentage categories, 500 base networks * 2000 simulations each
#
#
## taking the mean and sd of multiple networks in each percentage category
outbreak_variation_mean = np.mean(outbreak_variation_merged, axis = 1)

outbreak_variation_sd = np.std(outbreak_variation_merged, axis = 1)
#
#
## to check if it worked - all should return true
np.mean(outbreak_variation[1]) == outbreak_variation_mean[1]
np.std(outbreak_variation[3]) == outbreak_variation_sd[3]
#
#
## creating x - axis variable
x = [0,10,20,30,40,50,60,70,80,90]
# 
## plotting sd around the mean
plt.errorbar(x, outbreak_variation_mean*100, outbreak_variation_sd*100, linestyle='None', marker='.')
plt.xlabel('Percentage of missing nodes')
plt.ylabel('Change in outbreak size')
plt.title('Immunisation - Degree Centrality')
plt.show()


# also to be saved and reported - mean outbreak size of strategy without any removals
outbreak_list2 = []
for j in range(len(immunised_graph_list)): # for all percentage brackets
    mean_outbreak_list = []
    for g in immunised_graph_list[j]: # for all graphs in one percentage bracket
        inner_outbreaks = []
        while len(inner_outbreaks) < 500:  # 100 simulations for each graph
            sim = eon.Gillespie_SIR(g, 0.95, 1)
            inner_outbreaks.append(sim[1][0] - sim[1][-1]) # S(T) - S(0) (Outbreak Size) (L&F definition)
        mean_outbreak_list.append(inner_outbreaks)
    outbreak_list2.append(mean_outbreak_list)
    
np.mean(outbreak_list2)
np.std(outbreak_list2)

######################################################################################################

# TRADITIONAL ROBUSTNESS - code based on Martin and Niemeyer (2019)

# Base functions
def compare_centrality_dicts_correlation(d1, d2, scipy_correlation=kendalltau):
    if set(d1) != set(d2):
        nodes = sorted(set(d1).intersection(set(d2)))
    else:
        nodes = sorted(d1)

    v1 = np.round([d1[x] for x in nodes], 12)
    v2 = np.round([d2[x] for x in nodes], 12)

    return scipy_correlation(v1, v2).correlation


def robustness_calculator_builder(centrality_measure, 
                                  comparison_function=compare_centrality_dicts_correlation):
    @wraps(centrality_measure)
    def f(g0, g1):
        return compare_centrality_dicts_correlation(centrality_measure(g0), centrality_measure(g1))
    return f


def estimate_robustness(measured_network, error_mechanism, robustness_calculator, iterations=50, 
                        return_values=False):
    measured_robustness = np.array([robustness_calculator(measured_network, 
                                                          error_mechanism(measured_network))
                                    for _ in range(iterations)])
    vals = measured_robustness if return_values else None
    return namedtuple("robustness_estimate", "mean, sd values")(measured_robustness.mean(),
                                                                measured_robustness.std(), vals)

# a function for node removals - based on Martin's code
def remove_nodes_uniform(graph, alpha):
    nodes_to_keep = random.sample(graph.nodes, round(floor(graph.number_of_nodes())*(1 - alpha)))
    return graph.subgraph(nodes_to_keep)

# run "import random" again if it throws an error - fixes it
import random
# Loop over all graphs and sub-graphs
percent = [0,.10,.20,.30,.40,.50,.60,.70,.80,.90]
degree_robustness_calculator = robustness_calculator_builder(nx.degree_centrality)  

# for mean
robustness_list_mean = []
for j, p in zip(range(len(graph_list2)), percent): # for all percentage brackets
    inner_robustness_list = []
    for i in range(len(graph_list2[j])): # for all graphs in one percentage bracket
        rob = degree_robustness_calculator(graph_list2[j][i], remove_nodes_uniform(graph_list2[j][i], alpha = p))
        #m,s,v = estimate_robustness(graph_list2[j][i], partial(remove_nodes_uniform, alpha = p), 
        #                              degree_robustness_calculator,)
        inner_robustness_list.append(rob)
    robustness_list_mean.append(inner_robustness_list)


robustness_mean_plot = np.nanmean(robustness_list_mean, axis = 1)
robustness_sd_plot = np.nanstd(robustness_list_mean, axis = 1) 


# removing the first entries - not required as it calculates Tau between same entries
robustness_mean_plot = np.delete(robustness_mean_plot, 0)
robustness_sd_plot = np.delete(robustness_sd_plot, 0)

# Plotting it 
x = [10,20,30,40,50,60,70,80,90]
plt.errorbar(x, robustness_mean_plot*100, robustness_sd_plot*100, linestyle='None', marker='.')
plt.xlabel('Percentage of missing nodes')
plt.ylabel('Estimated Robustness (normalized)')
plt.title('Robustness - Degree Centrality')
plt.show()




# plot in the same figure with different colours
x1 = [0,10,20,30,40,50,60,70,80,90]
x2 = [10,20,30,40,50,60,70,80,90]
fig, ax = plt.subplots()
ax.errorbar(x1, outbreak_variation_mean*100, outbreak_variation_sd*100, linestyle='-', marker='.')
ax.errorbar(x2, robustness_mean_plot*100, robustness_sd_plot*100, linestyle='-', marker='.')

ax.set_xlabel('Percentage of missing nodes')
ax.set_ylabel('Change in outbreak size (% of pop.) \n vs \n Kendalls Tau (normalized)')
ax.set_title('Immunisation Robustness (blue) vs Traditional Robustness (orange) - Degree Centrality')

