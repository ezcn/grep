class Counter:
    counter = 0        #serve per incrementare i valori altrimenti prenderebve solo order(es.16)

    def __init__(self, init_value = 0):
        self.counter = init_value

    def inc(self, i=1):
        self.counter += i

    def dec(self, i=1):
        self.counter -= 1

def topological_sort(edges, source):
    order = Counter(len(edges))
    visited = [False] * len(edges)
    order_D_list = [0] * len(edges)

    recursive_topological_sort(edges, source, visited, order_D_list, order)     

    return order_D_list

def recursive_topological_sort(edges, v, visited, order_D_list, order):
    visited[v] = True
    for neighbour in edges[v]:                                 #per ogni vicino nel grafo partendo da v, se non visitato richiamare funzione
        if not visited[neighbour]:
            recursive_topological_sort(edges, neighbour, visited, order_D_list, order)
    
    order_D_list[v] = order.counter    
    order.dec()                        #decremento


def entrance(edges, v):
    child_one_parent_list = [True] * len(edges[v])

    for i, my_child in enumerate(edges[v]):
        for j, other_child_list in enumerate(edges):
            if my_child in other_child_list and v != j:
                child_one_parent_list[i] = False
                break

    return sum(child_one_parent_list) > 0   #se c'è almeno un figlio che appartiene ad un solo padre stampalo
                                            #altrimenti darà False

candidates_double_list = []

#insertentrance takes as input a vertex v, inserts it at the end of candidates and labels it as entrance
def insertentrance(v, candidates_double_list):
    candidates_double_list.append[(v,'entrance')]      

    return candidates_double_list['entrance']

#insertexit takes as input a vertex v, inserts it at the end of candidates and labels it as exit

def insertexit(v,candidates_double_list):
    candidates_double_list.append([v, 'exit'])

    return candidates_double_list['exit']

#head return the first candidate of a list
def head(candidates_double_list):    
    
    return candidates_double_list[0]

#return the last candidate of a list
def tail(candidates_double_list):   
    
    return candidates_double_list[-1]

#remove last candidate of a list
def delete_tail(candidates_double_list):
    del candidates_double_list[-1]         #si può anche togliere con candidate.list.pop

#return a follow candidate of v
def next(v,candidates_double_list):
    for i,(vertex,labels)in enumerate(candidates_double_list):
        if vertex == v:
            return candidates_double_list [i +1:]
            return[] 





'''
topological_sort(edges)
preventEnt = NULL
for v in topological_sort order
    alternativeEntrance[v] = NULL
    previousEntrance[v] = preventEnt
    if exit(v):
        insertexit(v)
    if entrance(v):
        insertentrance(v)
        preventEnt = v
while candidate is not empty   #lista vuota?
    if entrance(tail(candidate)):
        delete_tail(candidate)
    else reportsuperbubble(head(candidate),tail(candidate))

reportsuperbubble(start,exit)
if (start == NULL and exit == NULL and
    order_D_list[start] > order_D_list[exit]):
    delete_tail(candidate_list)
    return previousEntrance[exit]
'''    
# Python3 program to find minimum edge  
# between given two vertex of Graph 
import queue  #importa libreria
  
# function for finding minimum  
# no. of edge using BFS  
def minEdgeBFS(edges, u, n): #definisce minEdgeBFS dove edges sono archi, u radice ed n nodi
      
    # visited[n] for keeping track  
    # of visited node in BFS  
    visited = [False] * n  #inizializza con False
  
    ordered_nodes = []     #crea una lista vuota

    # Initialize distances as 0  
    distance = [0] * n 
  
    # queue to do BFS.  
    Q = queue.Queue() # crea una coda vuota
    #distance[u] = 0 # la distanza da se stesso è 0
  
    Q.put(u)  # Put the root-node in the queue
    visited[u] = True # Ho già visitato la radice (visto che parto da essa)
    while (not Q.empty()): # finché la coda non è vuota (finché ci sono vicini da visitare)
        x = Q.get()  # prende il nodo di cui visitare i vicini
          
        #print('x:', x)
        ordered_nodes.append(x)      #riempie la lista creata prima vuota, se i vicini del nodo stesso sono già stati visitati andare avanti
        for neighbour in edges[x]:   #per ogni vicino del nodo corrente se sono stati visitati i vicini, vai avanti
        	if visited[neighbour]:   #altrimenti calcolo la distanza
        		continue

        	print('neighbour:', neighbour)

        	# update distance for neighbour
        	distance[neighbour] = distance[x] + 1   #distanza del vicino è data da 0 più 1
        	Q.put(neighbour)
        	visited[neighbour] = True                #visita i vicini del nodo corrente, nodo corrente è presente nella coda
        	print('queue:', list(Q.queue))
        	print('visited:', visited)
        	print()
        
    return distance, ordered_nodes    #dà la distanza

# function for addition of edge  
def addEdge(edges, u, v): 
    edges[u].append(v)    # u --> v    #grafi diretti
    #edges[v].append(u)   # v --> u
  
# Driver  Code 
if __name__ == '__main__': 
  
  	# Num. of nodes
    n = 16
    
    # To store adjacency list of graph
    edges = []   #lista vuota
    for i in range(n):  # per ogni i nel range di n= nodi
    	edges.append([])  #riempi la lista con i nodi
 
    
	# Complete graph
    addEdge(edges, 1, 2)
    addEdge(edges, 2, 3)
    addEdge(edges, 1, 3)

    addEdge(edges, 3, 4)
    addEdge(edges, 3, 5)
    addEdge(edges, 3, 11)
    
    addEdge(edges, 5, 6)
    addEdge(edges, 5, 9)
    addEdge(edges, 9, 10)
    addEdge(edges, 6, 10)
    addEdge(edges, 10, 7)
    addEdge(edges, 6, 7)
    addEdge(edges, 7, 8)

    addEdge(edges, 4, 8)
    addEdge(edges, 11, 12)
    addEdge(edges, 12, 8)

    addEdge(edges, 8, 13)
    addEdge(edges, 8, 14)
    addEdge(edges, 13, 14)
    addEdge(edges, 13, 15)
    addEdge(edges, 15, 14)


    print('Edges')
    for i, edge in enumerate(edges):  #printa ogni nodo a quale altro nodo è legato, prendendo in considerazione l'indice
    	print(i, '-->', edge)
    print()


    print(topological_sort(edges, 1))
  
    u = 0
    distances_from_u, ordered_nodes = minEdgeBFS(edges, u, n)  #printa distanza
    print('Node and distance from root')
    for i, distance in enumerate(distances_from_u):  #lo stesso di sopra per la distanza
    	print(i, '-->', distance)
    print()

    dist_to_num_nodes = dict()    #dizionario che ha distanza dalla radice ordinata, se distanza unique già è nel dizionario
    for node, distance in enumerate(distances_from_u):  #enumerate: indice ed elemento, indice è uguale al nodo
    	if distance not in dist_to_num_nodes.keys():    #se la distanza non è tra le chiavi, metti 0, se già inizializzata aggiungi 1
    		dist_to_num_nodes[distance] = 0
    	dist_to_num_nodes[distance] += 1

    print('Distance from root --> Num. nodes')   #conto distanze 
    for k, v in dist_to_num_nodes.items():
    	print(k, '-->', v)

    start = True
    for node in ordered_nodes:        #per ogni distanza nella lista di nodi ordinati, entra usando la chiave distanza
   		key = distances_from_u[node]
   		if dist_to_num_nodes[key] == 1:  #se la distanza è univoca printa start altrimenti no
   			if start:
   				print(node, 'START')
   			else:
   				print(node, 'END')
   			start = not start
