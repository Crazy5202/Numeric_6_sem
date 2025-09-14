import os
import math
from copy import deepcopy

class ROTATE_SOLVER:
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

    def check_conditions(self) -> bool:
        """Проверить симметричность."""
        symm = True
        for i in range (self.n-1):
            for j in range(self.n-1):
                if (self.A[i][j]!=self.A[j][i]):
                    symm = False
            if not symm:
                print("Матрица не симметрична.")
                return False
        return True

    def stop_condition(self, current) -> float:
        """Посчитать корень суммы квадратов внедиагональных элементов."""
        sum = 0
        for i in range(self.n):
            for j in range(self.n):
                if (j>=i):
                    break
                sum += current[i][j]**2
        sum **= 0.5
        return sum

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
    
    def check_solution(self, vectors, values) -> bool:
        """Проверить решение на правильность."""
        for i in range(len(vectors)):
            left_multiply = [0.0]*self.n
            for j in range(self.n):
                for k in range(self.n):
                    left_multiply[j] += self.A[j][k]*vectors[i][k]
            
            right_multiply = [0.0]*self.n
            for j in range(self.n):
                right_multiply[j] = values[i]*vectors[i][j]

            for j in range (len(left_multiply)):
                if abs(left_multiply[j]-right_multiply[j]) > self.precision:
                    print("\nРЕШЕНИЕ НЕВЕРНО")
                    return False
        return True 

    def solve(self) -> tuple[list,list,int]:
        """Найти СЗ и СВ методом вращений."""
        if self.check_conditions() == False:
            return [], [], 0
        
        A_k = deepcopy(self.A)
        U = [[0.0]*self.n for _ in range (self.n)]
        for i in range(self.n):
            U[i][i] = 1
        counter = 0
        while (self.stop_condition(A_k)>self.precision):
            counter +=1
            max_val = A_k[0][1]
            max_ind = [0, 1]
            for i in range (self.n):
                for j in range(i):
                    if abs(A_k[i][j]) > abs(max_val):
                        max_val = A_k[i][j]
                        max_ind = [i,j]
            angle = 1/2*math.atan(2*A_k[max_ind[0]][max_ind[1]]/(A_k[max_ind[0]][max_ind[0]]-A_k[max_ind[1]][max_ind[1]]))
            U_k = [[0.0]*self.n for _ in range (self.n)]
            for i in range(self.n):
                U_k[i][i] = 1
            U_k[max_ind[0]][max_ind[0]] = math.cos(angle)
            U_k[max_ind[0]][max_ind[1]] = -math.sin(angle)
            U_k[max_ind[1]][max_ind[0]] = math.sin(angle)
            U_k[max_ind[1]][max_ind[1]] = math.cos(angle)
            U = self.multiply_matrix(U, U_k)
            A_k = self.multiply_matrix(self.transpose_matrix(U_k), A_k)
            A_k = self.multiply_matrix(A_k, U_k)
        values = []
        for i in range(self.n):
            values.append(A_k[i][i])
        vectors = []
        for i in range(self.n):
            vectors.append([U[j][i] for j in range(self.n)])
        if (self.check_solution(vectors, values) == False):
            return [], [], 0
        return values, vectors, counter
                

if __name__ == "__main__":
    matrix = ROTATE_SOLVER()

    round_num = matrix.read_from_file("matrix_4.txt")

    values, vectors, iters = matrix.solve()
    print("\nСобственные значения")
    for i in range(len(values)):
        print(f"СЗ_{i+1}: {round(values[i],round_num)}")

    print("\nСобственные векторы")
    for i in range(len(vectors)):
        print(f"СВ_{i+1}: ", end = "")
        for elem in vectors[i]:
            print(f"{round(elem,round_num)} ", end = "")
        print("")
        
    print(f"\nЧисло итераций: {iters}\n")