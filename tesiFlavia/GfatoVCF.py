import sys
sys.path.append('/home/flavia/Lab/odgi/lib')  
import odgi
g = odgi.graph()
g.load("/home/flavia/Lab/vgpop/lil.odgi")


start_node = g.get_handle(1)
g.get_id(start_node), g.get_sequence(start_node)

g_dfs = odgi.graph() # Array to keep track of visited nodes.

def create_edge_and_so_on(handle1, handle2, so_on_function, *args):
    handle1_id = g.get_id(handle1)
    handle2_id = g.get_id(handle2)
    if not g_dfs.has_node(handle2_id):
        so_on_function(args[0])

        if not g_dfs.has_node(handle1_id):
            g_dfs.create_handle(
                g.get_sequence(handle1),
                handle1_id
            )

        if not g_dfs.has_node(handle2_id):
            g_dfs.create_handle(
                g.get_sequence(handle2),
                handle2_id
            )
        g_dfs.create_edge(handle1, handle2)
        print('\tNew edge:', handle1_id, '-->', handle2_id)


def dfs(node_id):
    current_node = g.get_handle(node_id)
    sequence_node = g.get_sequence(current_node)
    print('current node_id:', node_id,'- sequence:', sequence_node)

    g.follow_edges(
        current_node,
        False,

        lambda neighbour:
            create_edge_and_so_on(
                current_node, neighbour, dfs, g.get_id(neighbour)
            )
    )

import queue


def calculate_distance(visited_node_id_set, prev_node_id, neighbour_id, Q, distances_dict):
    if neighbour_id not in visited_node_id_set:
        #print('neighbour_id:', neighbour_id)

        # update distance for neighbour
        distances_dict[neighbour_id] = distances_dict[prev_node_id] + 1
        Q.put(neighbour_id)
        visited_node_id_set.add(neighbour_id)
        #print('queue:', list(Q.queue))
        #print('visited_node_id_set:', visited_node_id_set)
        #print()

def bfs_distances(graph, starting_node_id):
    visited_node_id_set = set()
    ordered_node_id_list = []


    # lambdas don't permit assignment
    distances_dict = {}
    node_id_list = []
    graph.for_each_handle(lambda h: node_id_list.append(graph.get_id(h)))
    for node_id in node_id_list:
        distances_dict[node_id] = 0
    node_id_list.clear()


    Q = queue.Queue()
  
    Q.put(starting_node_id)
    visited_node_id_set.add(starting_node_id)
    while not Q.empty():
        current_node_id = Q.get()
        current_node = g.get_handle(current_node_id)

        #print('current_node_id:', current_node_id)
        ordered_node_id_list.append(current_node_id)
        
        graph.follow_edges(
            current_node,
            False,

            lambda neighbour:
                calculate_distance(
                    visited_node_id_set, current_node_id, graph.get_id(neighbour), Q, distances_dict
                )
        )

    return distances_dict, ordered_node_id_list

dfs(1)

def show_edge(a, b):
    print(g_dfs.get_id(a), "-->", g_dfs.get_id(b))

def display_node_edges(h):
    print("node", g_dfs.get_id(h))
    g_dfs.follow_edges(
        h, False,
        lambda n:
        show_edge(h, n))
    #g_dfs.follow_edges(
    #    h, True,
    #    lambda n:
    #    show_edge(n, h))
    
# displays all the edges twice, once for each of their ends
print('\nDFS graph')
g_dfs.for_each_handle(display_node_edges)
print()

distances_dict, ordered_node_id_list = bfs_distances(g_dfs, 1)

for node_id, distance in distances_dict.items():
    print(node_id, '- distance from root:', distance)

#print('ordered_node_id_lis:t', ordered_node_id_list)

dist_to_num_nodes = dict()    #dizionario che ha distanza dalla radice ordinata, se distanza unique già è nel dizionario
for node_id, distance in distances_dict.items():
    if distance not in dist_to_num_nodes.keys():    #se la distanza non è tra le chiavi, metti 0, se già inizializzata aggiungi 1
        dist_to_num_nodes[distance] = 0
    dist_to_num_nodes[distance] += 1

print('\nDistance from root --> Num. nodes')   #conto distanze 
for k, v in dist_to_num_nodes.items():
    print(k, '-->', v)

print('\nBubbles?')
start = True
for node_id in ordered_node_id_list:        #per ogni distanza nella lista di nodi ordinati, entra usando la chiave distanza
    key = distances_dict[node_id]
    if dist_to_num_nodes[key] == 1:  #se la distanza è univoca printa start altrimenti no
        if start:
            print(node_id, 'START') 
        else:
            print(node_id, 'END')
            print(get_sequences, 'seq')
            #print(node_id, 'START')
        start = not start
    if dist_to_num_nodes[key] != 1:
            print(node_id, 'Bolla')   
            print(sequence_node, "seq")

