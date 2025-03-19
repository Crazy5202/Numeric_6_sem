import os
import math
from copy import deepcopy

class QR_SOLVER:
    def __init__(self):
        self.A = []
        self.n = 0
        self.precision = 0

    def get_precision_num(self) -> int:
        """Получить число нулей погрешности."""
        prec_string = str(self.precision).split('.')
        if prec_string[0] != '0':
            if len(prec_string)==1:
                ret_val = int(prec_string[0].split('-')[1])
            else:
                ret_val = int(prec_string[1].split('-')[1])
            return ret_val
        count = 1
        for symb in prec_string[1]:
            if symb == '0':
                count += 1
            else:
                break
        return count

    def read_from_file(self, name) -> int:
        """Ввести матрицу из файла, возвращает число знаков для округления."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()
            self.n = len(lines)-1

            for i in range(self.n):
                line_stripped = lines[i].strip().split(' ')
                row = []
                for j in range (self.n):
                    row.append(float(line_stripped[j]))
                self.A.append(row)
            self.precision = float(lines[self.n])
        ret_val = self.get_precision_num()
        return ret_val

    def stop_condition(self, current) -> bool:
        """Посчитать корень суммы квадратов внедиагональных элементов."""
        for i in range(self.n):
            sum = 0
            for j in range(i+1, self.n):
                sum += current[j][i]**2
            sum **= 0.5
            if (sum > self.precision):
                return False
        return True

    def multiply_matrix(self, left, right) -> list:
        """Умножить матрицы."""
        multiplied = [[0.0] * len(right[0]) for _ in range(len(left))]
        for i in range (len(left)):
            for j in range(len(right[0])):
                for k in range(len(right)):
                    multiplied[i][j] += left[i][k]*right[k][j]
        return multiplied
    
    def transpose_matrix(self, matrix):
        transposed = [[0.0]*self.n for _ in range (self.n)]
        for i in range (self.n):
            for j in range(self.n):
                transposed[i][j] = matrix[j][i]
        return transposed

    def solve(self) -> tuple[list,int]:
        """Найти СЗ методом QR-разложения."""
        A_k = deepcopy(self.A)
        counter = 0
        while (not self.stop_condition(A_k)):
            counter += 1
            Q = [[0.0]*self.n for _ in range (self.n)]
            for i in range(self.n):
                Q[i][i] = 1
            for k in range (self.n):
                v = [0.0]*self.n
                v[k] = A_k[k][k] + A_k[k][k]/abs(A_k[k][k])*sum(A_k[i][k]**2 for i in range(k,self.n))**0.5
                for i in range (k+1, self.n):
                    v[i] = A_k[i][k]
                min_coeff = -2/sum(v[i]**2 for i in range (self.n))
                H = [[0.0]*self.n for _ in range (self.n)]
                for i in range (self.n):
                    for j in range (self.n):
                        H[i][j] = min_coeff*v[i]*v[j]
                    H[i][i] += 1
                A_k = self.multiply_matrix(H, A_k)
                Q = self.multiply_matrix(Q, H)
            A_k = self.multiply_matrix(A_k, Q)
        values = [A_k[i][i] for i in range(self.n)]
        return values, counter          

if __name__ == "__main__":
    matrix = QR_SOLVER()

    round_num = matrix.read_from_file("matrix_5.txt")

    values, iters = matrix.solve()
    print("\nСобственные значения")
    for i in range(len(values)):
        print(f"СЗ_{i+1}: {round(values[i],round_num)}")

    print(f"\nЧисло итераций: {iters}\n")