import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment
from itertools import product
from grid import Grid


class Solver:
    """
    A solver class. 

    Attributes: 
    -----------
    grid: Grid
        The grid
    pairs: list[tuple[tuple[int]]]
        A list of pairs, each being a tuple ((i1, j1), (i2, j2))
    """

    def __init__(self, grid):
        """
        Initializes the solver.

        Parameters: 
        -----------
        grid: Grid
            The grid
        """
        self.grid = grid
        self.pairs = list()

    def score(self):
        """
        Computes the score of the list of pairs in self.pairs
        """
        total_score = 0
        used_cells = set()

        # Cellules prises dans une paire
        for (i1, j1), (i2, j2) in self.pairs:
            total_score += self.grid.cost(((i1, j1), (i2, j2)))
            used_cells.add((i1, j1))
            used_cells.add((i2, j2))

        # Cellules non prises dans une paire (hors cellules noires)
        for i in range(self.grid.n):
            for j in range(self.grid.m):
                if (i, j) not in used_cells and not self.grid.is_forbidden(i, j):
                    total_score += self.grid.value[i][j]

        return total_score

class SolverEmpty(Solver):
    
    def run(self):
        pass

class SolverGreedy(Solver): # Complexité totale : O(nm log(nm))

    def run(self):
        self.pairs = []
        taken_cells = set()

        # Complexité : O(nm)
        valid_pairs = self.grid.all_pairs()

        # Tri des paires par coût croissant (Complexité : O(nm log(nm)))
        valid_pairs.sort(key=lambda pair: self.grid.cost(pair))

        # Complexité : O(nm)
        for pair in valid_pairs:
            (i1, j1), (i2, j2) = pair
            if (i1, j1) in taken_cells or (i2, j2) in taken_cells:
                continue
            self.pairs.append(pair)
            taken_cells.add((i1, j1))
            taken_cells.add((i2, j2))

        return self.pairs

class SolverBipartiteMatching(Solver): # Complexité : O(nm sqrt(nm))
    def run(self):
        G = nx.Graph()
        pairs = self.grid.all_pairs()
        left_nodes = set()
        right_nodes = set()
        
        # Graphe biparti (Complexité : O(nm))
        for (i1, j1), (i2, j2) in pairs:
            if (i1 + j1) % 2 == 0: # Pair
                left_nodes.add((i1, j1))
                right_nodes.add((i2, j2))
            else: # Impair
                left_nodes.add((i2, j2))
                right_nodes.add((i1, j1))
            G.add_edge((i1, j1), (i2, j2))
        
        # Matching maximal (Returns the maximum cardinality matching in the given bipartite graph) (Complexité : O(nm sqrt(nm)))
        matching = nx.bipartite.hopcroft_karp_matching(G, top_nodes=left_nodes)

        # Filtrage (Complexité : O(nm))
        self.pairs = []
        taken_nodes = set()
        for u, v in matching.items():
            if u in left_nodes and v in right_nodes and u not in taken_nodes and v not in taken_nodes:
                self.pairs.append((u, v))
                taken_nodes.add(u)
                taken_nodes.add(v)
        
        return self.pairs

class SolverBipartiteMatching_bis(Solver):
    def run(self):
        pairs = self.grid.all_pairs()
        left_nodes = set()
        right_nodes = set()
        adj = {}

        # Construire le graphe biparti sous forme de dictionnaire d'adjacence (O(nm))
        for (i1, j1), (i2, j2) in pairs:
            if (i1 + j1) % 2 == 0:  # Pair → Groupe gauche
                left_nodes.add((i1, j1))
                right_nodes.add((i2, j2))
                adj.setdefault((i1, j1), []).append((i2, j2))
                adj.setdefault((i2, j2), []).append((i1, j1))
            else:  # Impair → Groupe gauche
                left_nodes.add((i2, j2))
                right_nodes.add((i1, j1))
                adj.setdefault((i2, j2), []).append((i1, j1))
                adj.setdefault((i1, j1), []).append((i2, j2))

        def hopcroft_karp(left_nodes, right_nodes, adj):
            """
            Implémente l'algorithme de Hopcroft-Karp pour trouver un appariement maximum dans un graphe biparti.
            Complexité : O(nm sqrt(nm))
            
            - left_nodes : ensemble des nœuds de gauche
            - right_nodes : ensemble des nœuds de droite
            - adj : dictionnaire d'adjacence
            """
            matching = {}
            distance = {}

            def bfs():
                """
                Phase de BFS pour trouver un chemin d'augmentation.
                Retourne True si un chemin d'augmentation est trouvé.
                """
                queue = []
                for u in left_nodes:
                    if u not in matching:  # Sommets non appariés
                        distance[u] = 0
                        queue.append(u)
                    else:
                        distance[u] = float('inf')

                distance[None] = float('inf')

                for u in queue:
                    if distance[u] < distance[None]:  # On continue jusqu'à trouver un chemin d'augmentation
                        for v in adj.get(u, []):
                            match_v = matching.get(v)
                            if distance.get(match_v, float('inf')) == float('inf'):
                                distance[match_v] = distance[u] + 1
                                queue.append(match_v)

                return distance[None] != float('inf')

            def dfs(u):
                """
                Phase de DFS pour améliorer l'appariement via les chemins d'augmentation trouvés par BFS.
                """
                if u is not None:
                    for v in adj.get(u, []):
                        match_v = matching.get(v)
                        if distance.get(match_v, float('inf')) == distance[u] + 1 and dfs(match_v):
                            matching[v] = u
                            matching[u] = v
                            return True
                    distance[u] = float('inf')
                    return False
                return True

            # Boucle principale de l'algorithme
            while bfs():
                for u in left_nodes:
                    if u not in matching:
                        dfs(u)

            return matching

        # Appliquer l'algorithme de Hopcroft-Karp
        matching = hopcroft_karp(left_nodes, right_nodes, adj)

        # Convertir le dictionnaire de matching en liste de paires
        self.pairs = [(u, v) for u, v in matching.items() if u in left_nodes and v in right_nodes]

        return self.pairs

# Question 4 :  
# Pour le graphe suivant :
# 2 3 
# 0 0 0 
# 0 0 0 
# 9 9 8
# 8 0 0
# Greedy solver ne renvoie pas la réponse optimale qui est de score 2 mais une de score 16
# On pourrait utiliser une méthode de Brute Force dont la complexité serait en O(2^(nm))

class SolverHungarian(Solver): # Complexité totale : O((nm)^3)
    def run(self):
        # Générer toutes les paires possibles pour une case (i, j)  (Complexité : O(nm))   
        def all_pairs_bis(i, j, n, m):
            return [[(i, j), (k, l)] for k in range(n) for l in range(m)]
        
        # Générer toutes les paires possibles pour toutes les cases (i, j) (Complexité : O(n^2 m^2))
        def generate_pairs_matrix(n, m):
            pairs_matrix = [all_pairs_bis(i, j, n, m) for i in range(n) for j in range(m)]
            return np.array(pairs_matrix, dtype=object)
        
        pairs_matrix = generate_pairs_matrix(self.grid.n, self.grid.m)

        # Convertir all_pairs en un ensemble de tuples hachables (Complexité : O((nm)^2))
        all_valid_pairs = set(map(tuple, self.grid.all_pairs()))

        # Générer la matrice de coût (Complexité : O((nm)^2))
        cost_matrix = np.array([
            [
                self.grid.cost(pair) if tuple(map(tuple, pair)) in all_valid_pairs
                else self.grid.value[pair[0][0]][pair[0][1]] + self.grid.value[pair[1][0]][pair[1][1]]
                for pair in row
            ] 
            for row in pairs_matrix
        ], dtype=float)

        for k in range(self.grid.n * self.grid.m): # Complexité : O((nm)^2)
            cost_matrix[k][k] = np.inf

        # Appliquer l'algorithme hongrois pour le problème d'appariement (Complexité : O((nm)^3))
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        # Lister les cellules de la grille
        cells = [(i, j) for i in range(self.grid.n) for j in range(self.grid.m)]

        # Associer les indices retournés par l’algorithme aux cellules de la grille (Complexité : O((nm)^2))
        self.pairs = []
        used_cells = set()

        for r, c in zip(row_ind, col_ind):
            if r < len(cells) and c < len(cells):  # Vérifier que les indices sont valides
                cell1 = cells[r]
                cell2 = cells[c]

                # Vérifier si la paire est valide et si les cellules ne sont pas déjà prises
                if (
                    (cell1, cell2) in all_valid_pairs and
                    cell1 not in used_cells and
                    cell2 not in used_cells
                ):
                    self.pairs.append((cell1, cell2))
                    used_cells.add(cell1)
                    used_cells.add(cell2)

        return self.pairs