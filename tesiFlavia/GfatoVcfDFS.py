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
    n = 10
    
    # To store adjacency list of graph
    edges = []   #lista vuota
    for i in range(n):  # per ogni i nel range di n= nodi
    	edges.append([])  #riempi la lista con i nodi
 
    '''
	# Complete graph
    addEdge(edges, 0, 1)
    addEdge(edges, 1, 2)
    addEdge(edges, 2, 3)
    addEdge(edges, 0, 5)
    addEdge(edges, 5, 4)
    addEdge(edges, 4, 3)
    addEdge(edges, 3, 6)
    addEdge(edges, 6, 7)
    addEdge(edges, 6, 8)
    addEdge(edges, 7, 9)
    addEdge(edges, 8, 9)
    '''

    # DFS tree
    addEdge(edges, 0, 1)
    addEdge(edges, 1, 2)
    addEdge(edges, 2, 3)
    addEdge(edges, 0, 5)
    addEdge(edges, 5, 4)
    #addEdge(edges, 4, 3)
    addEdge(edges, 3, 6)
    addEdge(edges, 6, 7)
    addEdge(edges, 6, 8)
    addEdge(edges, 7, 9)
    #addEdge(edges, 8, 9)

    print('Edges')
    for i, edge in enumerate(edges):  #printa ogni nodo a quale altro nodo è legato, prendendo in considerazione l'indice
    	print(i, '-->', edge)
    print()
  
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

'''
non usato
visited = [] # Array to keep track of visited nodes.
def dfs(visited, edges, node):
    if node not in visited:          #se il nodo non è già visitato printa il nodo e aggiungilo alla lista, per ogni vicino del nodo
        print(node)
        visited.append(node)
        for neighbour in edges[node]:
            dfs(visited, edges, neighbour)

#dfs(visited, edges, 0)
'''