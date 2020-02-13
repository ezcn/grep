# Create a graph
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms.traversal.depth_first_search import dfs_tree
from networkx.algorithms.traversal.breadth_first_search import bfs_tree

#graph = Graph()
graph = nx.Graph() #grafo vuoto diretto
# Paper-ino
graph.add_edge(0, 1)
graph.add_edge(0, 5)
graph.add_edge(1, 2)
graph.add_edge(1, 6)
graph.add_edge(2, 3)
graph.add_edge(2, 4)
graph.add_edge(3, 4)
graph.add_edge(4, 5)
graph.add_edge(4, 1)
graph.add_edge(5, 6)
graph.add_edge(6, 7)

    # graph.add_edge(0, 1)
    # graph.add_edge(0, 2)
    # graph.add_edge(2, 1)
    # graph.add_edge(1, 3)
    # graph.add_edge(2, 3)
    # graph.add_edge(3, 0)




nx.draw(graph, node_size=500) #visualizza grafo
plt.show()