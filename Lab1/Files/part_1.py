import os
from copy import deepcopy
from pprint import pprint

class LU_SOLVER:
    def __init__(self):
        self.A = []
        self.b = []
        self.n = 0
        self.LU = []

    def read_from_file(self, name):
        """Ввести матрицу из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()
            self.n = len(lines)

            for i in range(self.n):
                line_stripped = lines[i].strip().split(' ')
                row = []
                for j in range (self.n):
                    row.append(float(line_stripped[j]))
                self.A.append(row)
                self.b.append(float(line_stripped[self.n]))

    def decompose(self) -> list:
        """Провести LU-разложение матрицы."""
        LU = deepcopy(self.A)

        for k in range(self.n):
            for i in range(k+1, self.n):
                LU[i][k] /= LU[k][k]
                for j in range(k+1, self.n):
                    LU[i][j] -= LU[i][k]*LU[k][j]

        self.LU = LU
        return LU

    def check_solution(self, right, checking) -> bool:
        """Проверить решение на правильность."""
        calc_b = []
        for i in range(self.n):
            new_val = 0
            for j in range (self.n):
                new_val += self.A[i][j]*checking[j]
            calc_b.append(new_val)
            
            if abs(calc_b[i]-right[i])>1e-3:
                print("Решение неверно")
                return False
        return True 

    def solve_func(self, right) -> list:
        """Решить СЛАУ через LU-разложение для заданной правой части."""
        if (self.LU==[]):
            self.decompose()
        z = []
        z.append(right[0])
        for i in range(1,self.n):
            new_val = right[i]
            for j in range (i):
                new_val -= z[j]*self.LU[i][j]
            z.append(new_val)
        x = deepcopy(z)
        for i in range(self.n-1, -1, -1):
            for j in range (i+1, self.n):
                x[i] -= self.LU[i][j]*x[j]
            x[i] /= self.LU[i][i]

        if (self.check_solution(right, x) == False):
            return []

        return x
    
    def solve(self) -> list:
        """Решить СЛАУ через LU-разложение для заданной матрицы."""
        return self.solve_func(self.b)

    def calc_determinant(self) -> float:
        """Найти детерминант матрицы из LU-разложения."""
        if self.LU == []:
            self.decompose()

        ans = self.LU[0][0]
        for i in range(1,self.n):
            ans *= self.LU[i][i]

        return ans
    
    def calc_inverse(self) -> tuple[list, list]:
        """Рассчитать обратную матрицу."""
        ans_inv = []
        for i in range(self.n):
            col = [0.0]*self.n
            col[i] = 1.0
            ans_inv.append(self.solve_func(col))
        ans = [[0.0]*self.n for _ in range (self.n)]
        for i in range (self.n):
            for j in range (self.n):
                ans[i][j] = ans_inv[j][i]
        multiplied = [[0.0]*self.n for _ in range (self.n)]
        for i in range (self.n):
            for j in range(self.n):
                for k in range(self.n):
                    multiplied[i][j] += self.A[i][k]*ans[k][j]
        return ans_inv, multiplied
        

if __name__ == "__main__":
    matrix = LU_SOLVER()

    matrix.read_from_file("matrix_1.txt")

    print("\nМатрица LU:")
    for row in matrix.decompose():
        for elem in row:
            if (elem>0 and elem < 10):
                print(" ", end="")
            print(f"{round(elem,1)} ", end = "")
        print("")

    print("\nРешение СЛАУ:")
    for elem in matrix.solve():
        print(round(elem,3), end = " ")

    print(f"\n\nОпределитель матрицы: {round(matrix.calc_determinant(),3)}")

    inverted, check = matrix.calc_inverse()
    print("\nОбратная матрица:")
    for row in inverted:
        for elem in row:
            print(f"{round(elem,3)} ", end = "")
        print("")
    
    print("Произведение на обратную матрицу:")
    for row in check:
        for elem in row:
            print(f"{round(elem,3)} ", end = "")
        print("")