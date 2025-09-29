import matplotlib.pyplot as plt
import numpy as np
import os

"""
This is the grid module. It contains the Grid class and its associated methods.
"""

class Grid():
    """
    A class representing the grid. 

    Attributes: 
    -----------
    n: int
        Number of lines in the grid
    m: int
        Number of columns in the grid
    color: list[list[int]]
        The color of each grid cell: value[i][j] is the value in the cell (i, j), i.e., in the i-th line and j-th column. 
        Note: lines are numbered 0..n-1 and columns are numbered 0..m-1.
    value: list[list[int]]
        The value of each grid cell: value[i][j] is the value in the cell (i, j), i.e., in the i-th line and j-th column. 
        Note: lines are numbered 0..n-1 and columns are numbered 0..m-1.
    colors_list: list[char]
        The mapping between the value of self.color[i][j] and the corresponding color
    """
    

    def __init__(self, n, m, color=[], value=[]):
        """
        Initializes the grid.

        Parameters: 
        -----------
        n: int
            Number of lines in the grid
        m: int
            Number of columns in the grid
        color: list[list[int]]
            The grid cells colors. Default is empty (then the grid is created with each cell having color 0, i.e., white).
        value: list[list[int]]
            The grid cells values. Default is empty (then the grid is created with each cell having value 1).
        
        The object created has an attribute colors_list: list[char], which is the mapping between the value of self.color[i][j] and the corresponding color
        """
        self.n = n
        self.m = m
        if not color:
            color = [[0 for j in range(m)] for i in range(n)]            
        self.color = color
        if not value:
            value = [[1 for j in range(m)] for i in range(n)]            
        self.value = value
        self.colors_list = ['w', 'r', 'b', 'g', 'k']

    def __str__(self): 
        """
        Prints the grid as text.
        """
        output = f"The grid is {self.n} x {self.m}. It has the following colors:\n"
        for i in range(self.n): 
            output += f"{[self.colors_list[self.color[i][j]] for j in range(self.m)]}\n"
        output += f"and the following values:\n"
        for i in range(self.n): 
            output += f"{self.value[i]}\n"
        return output

    def __repr__(self): 
        """
        Returns a representation of the grid with number of rows and columns.
        """
        return f"<grid.Grid: n={self.n}, m={self.m}>"

    def plot(self): 
        """
        Plots a visual representation of the grid.
        """
        color_map = {0: 'white', 1: 'red', 2: 'blue', 3: 'green', 4: 'black'}
        n, m = self.n, self.m
        fig, ax = plt.subplots(figsize=(m, n))

        for i in range(n):
            for j in range(m):
                cell_color = color_map[self.color[i][j]]
                rect = plt.Rectangle((j, n - 1 - i), 1, 1, facecolor=cell_color, edgecolor='black')
                ax.add_patch(rect)
                ax.text(j + 0.5, n - 1 - i + 0.5, str(self.value[i][j]), 
                        color='black' if self.color[i][j] != 4 else 'white',
                        fontsize=12, ha='center', va='center')
        
        ax.set_xlim(0, m)
        ax.set_ylim(0, n)
        ax.set_xticks(np.arange(0, m, 1))
        ax.set_yticks(np.arange(0, n, 1))
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(True, color='black', linewidth=1)
        
        plt.show()

    def is_forbidden(self, i, j):
        """
        Returns True is the cell (i, j) is black and False otherwise
        """
        return self.color[i][j] == 4

    def cost(self, pair):
        """
        Returns the cost of a pair
 
        Parameters: 
        -----------
        pair: tuple[tuple[int]]
            A pair in the format ((i1, j1), (i2, j2))

        Output: 
        -----------
        cost: int
            the cost of the pair defined as the absolute value of the difference between their values
        """
        (i1, j1), (i2, j2) = pair
        return abs(self.value[i1][j1] - self.value[i2][j2])

    def are_compatible(self, cell1, cell2):
        """
        Returns True if cell1 and cell2 can be paired together.
        """
        i1, j1 = cell1
        i2, j2 = cell2
        c1, c2 = self.color[i1][j1], self.color[i2][j2]

        compatible_colors = {
            0: {1, 2, 3, 0}, # Blanc : Blanc, Rouge, Bleu, Vert et Noir
            1: {1, 2, 0}, # Rouge : Blanc, Rouge et Bleu
            2: {1, 2, 0}, # Bleu : Blanc, Rouge et Bleu
            3: {3, 0}, # Vert : Blanc ou Vert
            4: {} # Noir :
        }

        return c2 in compatible_colors[c1]

    def all_pairs(self):
        """
        Returns a list of all pairs of cells that can be taken together. 

        Outputs a list of tuples of tuples [(c1, c2), (c1', c2'), ...] where each cell c1 etc. is itself a tuple (i, j)
        """
        pairs = []
        for i in range(self.n):
            for j in range(self.m):
                if self.is_forbidden(i, j):
                    continue  
                
                # Test de la case de droite
                if j + 1 < self.m and not self.is_forbidden(i, j + 1):
                    if self.are_compatible((i, j), (i, j + 1)):
                        pairs.append(((i, j), (i, j + 1)))
                
                # Test de la case du bas
                if i + 1 < self.n and not self.is_forbidden(i + 1, j):
                    if self.are_compatible((i, j), (i + 1, j)):
                        pairs.append(((i, j), (i + 1, j)))

        return pairs

    @classmethod
    def grid_from_file(cls, file_name, read_values=False): 
        """
        Creates a grid object from class Grid, initialized with the information from the file file_name.
        
        Parameters: 
        -----------
        file_name: str
            Name of the file to load. The file must be of the format: 
            - first line contains "n m" 
            - next n lines contain m integers that represent the colors of the corresponding cell
            - next n lines [optional] contain m integers that represent the values of the corresponding cell
        read_values: bool
            Indicates whether to read values after having read the colors. Requires that the file has 2n+1 lines

        Output: 
        -------
        grid: Grid
            The grid
        """
        with open(file_name, "r") as file:
            n, m = map(int, file.readline().split())
            color = [[] for i_line in range(n)]
            for i_line in range(n):
                line_color = list(map(int, file.readline().split()))
                if len(line_color) != m: 
                    raise Exception("Format incorrect")
                for j in range(m):
                    if line_color[j] not in range(5):
                        raise Exception("Invalid color")
                color[i_line] = line_color

            if read_values:
                value = [[] for i_line in range(n)]
                for i_line in range(n):
                    line_value = list(map(int, file.readline().split()))
                    if len(line_value) != m: 
                        raise Exception("Format incorrect")
                    value[i_line] = line_value
            else:
                value = []

            grid = Grid(n, m, color, value)
        return grid