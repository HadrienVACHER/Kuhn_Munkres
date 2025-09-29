from grid import Grid
from solver import *
import os

data_path = "../input/"

file_name = data_path + "grid04.in"
grid = Grid.grid_from_file(file_name, read_values=True)

solver = SolverEmpty(grid)
solver.run()
print("The final score of SolverEmpty is:", solver.score())

solver = SolverGreedy(grid)
solver.run()
print("The final score of SolverGreedy is:", solver.score())

solver = SolverBipartiteMatching(grid)
solver.run()
print("The final score of SolverBipartite is:", solver.score())

solver = SolverBipartiteMatching_bis(grid)
solver.run()
print("The final score of SolverBipartite_bis is:", solver.score())

solver = SolverHungarian(grid)
solver.run()
print("The final score of SolverHungarian is:", solver.score())

grid.plot()


grid_files = {
    "grid00.in": 12,
    "grid01.in": 8,
    "grid02.in": 1,
    "grid03.in": 2,
    "grid04.in": 4,
    "grid05.in": 35,
    "grid11.in": 26,
    "grid12.in": 19,
    "grid13.in": 22,
    "grid14.in": 27,
    "grid15.in": 21,
    "grid16.in": 28,
    "grid17.in": 256,
    "grid18.in": 259,
    "grid19.in": 248
}

solvers = {
    "Empty": SolverEmpty,
    "Greedy": SolverGreedy,
    "BipartiteMatching": SolverBipartiteMatching,
    "BipartiteMatching_bis": SolverBipartiteMatching_bis,
    "Hungarian": SolverHungarian
}

# Tableau des résultats
print("{:<15} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}".format("Grid", "Empty", "Greedy", "Bipartite", "Bipartite_bis", "Hungarian", "Score Optimal"))
print("-" * 80)

for file_name, expected_score in grid_files.items():
    grid = Grid.grid_from_file(os.path.join(data_path, file_name), read_values=True)
    scores = {}
    
    for solver_name, SolverClass in solvers.items():
        solver = SolverClass(grid)
        solver.run()
        scores[solver_name] = solver.score()
    
    print("{:<15} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10}".format(
        file_name, scores["Empty"], scores["Greedy"], scores["BipartiteMatching"], scores["BipartiteMatching_bis"], scores["Hungarian"], expected_score
    ))


# SolverHungarian ne renvoie pas la réponse optimale mais il reste quand même pas si mal, du moins sur les grandes grilles