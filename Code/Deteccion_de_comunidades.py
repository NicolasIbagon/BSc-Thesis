
#El preprocseamiento y la ejecucion de parte de los algoritmos de deteccion de comunidades se realiza en este programa

import matplotlib.pyplot as plt
#import numpy as np
import pandas as pd
import networkx as nx
import csv
from sys import stdin
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community import k_clique_communities
from networkx.algorithms.community import asyn_lpa_communities
from networkx.algorithms.community import modularity

#import markov_clustering as mc
#import igraph
#conda install -c conda-forge python-levenshtein
#from cdlib import algorithms
#import cairocffi
#import igraph

filename = "C:/Users/HOME/Desktop/Universidad Javeriana Cali/Tesis/pnas.2025581118.sd02.csv"
#D:/PROGRAMAS/Dropbox/Tesis Ingeniería Cancer de Colon/pnas.2025581118.sd02.csv

data = pd.read_csv(filename)
Graphtype = nx.Graph()
G = nx.from_pandas_edgelist(data, 
                            source="proteinA_entrezid",
                            target="proteinB_entrezid",
                            edge_attr='databases', create_using=Graphtype)

G.remove_edges_from(list(nx.selfloop_edges(G)))
G2 = nx.to_undirected(G)


Gcc = sorted(nx.connected_components(G2), key=len, reverse=True)

G0 = G.subgraph(Gcc[0]) #Hay 18446 nodos en G0
#---------------------------------------------------------------------------------------


a = []
b = []



ent = stdin.readline().strip() #------Leer nodos
while ent != "":
	a.append(int(ent))
	ent = stdin.readline().strip()
print("a ready")

ent = stdin.readline().strip() #------Leer valores ponderados
while ent != "":
	b.append(float(ent))
	ent = stdin.readline().strip()
print("b ready")

for node in range(len(a)): #----------Asignar valores ponderados a los nodos
	if a[node] in G0:
		G0.nodes[a[node]]["pond"] = b[node]


"""
for i in range(len(array_nodes_0)):
	G0.remove_node(array_nodes_0[i])
"""



for edge in G0.edges: #---------------Asignacion de los pesos de las conexiones a los ejes 
	G0[edge[0]][edge[1]]["weight"] = (G.nodes[edge[0]]["pond"]) + (G.nodes[edge[1]]["pond"])



#Se borran los nodos con ponderado 0

G0_unfrozen= nx.Graph(G0)
array_edges_0 = []
for edge in G0_unfrozen.edges:
	if G0_unfrozen[edge[0]][edge[1]]["weight"] == 0:
		array_edges_0.append((edge[0], edge[1]))

G0_unfrozen.remove_edges_from(array_edges_0)
isolates = nx.isolates(G0_unfrozen)

G0_unfrozen.remove_nodes_from(list(isolates))

for edge in sorted(G0_unfrozen.edges, key=lambda x:G0_unfrozen[x[0]][x[1]]["weight"],  reverse=True):              #--------- imprimir pesos de los edges
	print(str(edge[0]), str(edge[1]), G0[edge[0]][edge[1]]["weight"])




#Algoritmos de detección de comunidades


dic = {}

"""
Girvan newman

first_communities = nx.algorithms.community.girvan_newman(G0_unfrozen)
node_groups = []
for com in next(first_communities):
  node_groups.append(com)

print(node_groups)

"""

"""
Asyn fluid
"""
second_communities = nx.algorithms.community.asyn_fluidc(G0, 100)
for idx, com in enumerate(second_communities):
    dic[str(idx)] = list(com)

"""
asyn_lpa_communities


third_communities = nx.algorithms.community.asyn_lpa_communities(G0_unfrozen, "weight")

for idx, com in enumerate(third_communities):
    dic[str(idx)] = list(com)
"""


"""
Greedy modularity communities



four_communities = nx.algorithms.community.greedy_modularity_communities(G0_unfrozen, "weight")
for idx, com in enumerate(four_communities):
    dic[str(idx)] = list(com)
"""





#Escribir resultador algoritmos comunidades en csv


file_name = 'asyn_fluidc.csv'
header = ['Id', 'modularity_class']
with open(file_name,'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)

    for key in dic:
        for i in range(len(dic[key])):
            writer.writerow([dic[key][i],key])

"""
ponderados = {}
for node in range(len(a)):
	ponderados[a[node]] = b[node]

print(ponderados[9796], ponderados[7918], ponderados[4609])"""


#for i in sorted(G2.degree, key=lambda x: x[1], reverse=True): #--------- Imprimir los ponderados en orden
#	print(i[0],i[1])


"""
centrality = nx.betweenness_centrality(G2, k=50, endpoints=True) #------- Betweeness centrality 

sorted_values = sorted(centrality.values()) # Sort the values
sorted_dict = {}

for i in sorted_values:
    for k in centrality.keys():
        if centrality[k] == i:
            sorted_dict[k] = centrality[k]
            break

print(sorted_dict)
"""




#for edge in G0.edges:              #--------- imprimir pesos de los edges
#	print(edge, G0[edge[0]][edge[1]]["weight"])

#G0[10485][100631383]["weight"] = 5 --------- asignar peso a un eje ya creado
#print(G0.edges)                    --------- lista de ejes
#print(G0[10485][100631383]["weight"]) ------ imprimir peso de eje
#G.nodes[1]["pond"] = 714  ------------------ establecer valor ponderado de nodo

#nx.write_edgelist(G2, "D:/USER/Trabajos/Trabajos/Trabajo de grado/Python/pnas_G_original.csv",data=False)
#nx.write_edgelist(G0, "D:/USER/Trabajos/Trabajos/Trabajo de grado/Python/pnas_GO.csv",data=False)
