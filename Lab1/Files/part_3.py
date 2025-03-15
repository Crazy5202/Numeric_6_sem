import os

class Iter:
    def __init__(self):
        self.A = []
        self.b = []
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
        
    def alpha_metric(self) -> float:
        """Посчитать норму матрицы альфа."""
        sum = 0
        for i in range(self.n):
            for j in range(self.n):
                if (i==j):
                    continue
                sum += (self.A[i][j]/self.A[i][i])**2
        sum **= 0.5
        return sum

    def error_metric(self, vec1, vec2) -> float:
        """Посчитать норму разности вектора решений текущей и прошлой итерации."""
        sum = 0
        for i in range (len(vec1)):
            sum += (vec1[i]-vec2[i])**2
        sum **= 0.5
        return sum

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
                self.b.append(float(line_stripped[self.n]))
            self.precision = float(lines[self.n])
        ret_val = self.get_precision_num()
        self.precision *= (1/self.alpha_metric()-1)
        return ret_val
        

    def check_conditions(self) -> bool:
        """Проверить сходимость."""
        for i in range (self.n):
            not_diag = 0
            for j in range(self.n):
                if (i==j):
                    continue
                not_diag += abs(self.A[i][j])
            if abs(self.A[i][i]) <= not_diag:
                print("Не удовлетворяет условию сходимости.")
                return False
        print("Удовлетворяет условию сходимости.")
        return True
    
    def check_solution(self, x) -> bool:
        """Проверить решение на правильность."""
        calc_b = []
        for i in range(self.n):
            new_val = 0
            for j in range (self.n):
                new_val += self.A[i][j]*x[j]
            calc_b.append(new_val)
            
            if abs(calc_b[i]-self.b[i])>self.precision/(1/self.alpha_metric()-1):
                print("Решение неверно")
                return False
        print("Решение верно")
        return True 

    def stop_iter(self, vec1, vec2) -> bool:
        if self.error_metric(vec1, vec2) <= self.precision:
            return True
        return False

    def simple_iter(self) -> tuple[list,int]:
        """Решить методом простых итераций."""
        if self.check_conditions() == False:
            return [], 0
        x = [self.b[i]/self.A[i][i] for i in range (self.n)]
        count = 0
        while True:
            count +=1
            new_x = []
            for i in range (self.n):
                new_val = self.b[i]/self.A[i][i]
                for j in range (self.n):
                    if i==j: continue
                    new_val = new_val - x[j]*self.A[i][j]/self.A[i][i]
                new_x.append(new_val)
            if self.stop_iter(new_x, x):
                if self.check_solution(new_x) == False:
                    return [], 0
                return new_x, count
            x = new_x
        
    
    def zeydel(self) -> tuple[list, int]:
        """Решить методом Зейделя."""
        if self.check_conditions() == False:
            return [], 0
        x = [self.b[i]/self.A[i][i] for i in range (self.n)]
        count = 0
        while True:
            count += 1
            new_x = x.copy()
            for i in range (self.n):
                new_x[i] = self.b[i]/self.A[i][i]
                for j in range (self.n):
                    if i==j: continue
                    new_x[i] -= new_x[j]*self.A[i][j]/self.A[i][i]
            if self.stop_iter(new_x, x):
                if self.check_solution(new_x) == False:
                    return [], 0
                return new_x, count
            x = new_x

if __name__ == "__main__":
    matrix = Iter()
    round_num = matrix.read_from_file("matrix_3.txt")
    matrix.get_precision_num()
    result, iters = matrix.simple_iter()
    print("\nМЕТОД ПРОСТЫХ ИТЕРАЦИЙ\n")
    for i in range (len(result)):
        print(f"x_{i+1} = {round(result[i],round_num)}\n")
    print(f"Число итераций: {iters}")
    result, iters = matrix.zeydel()
    print("\nМЕТОД ЗЕЙДЕЛЯ\n")
    for i in range (len(result)):
        print(f"x_{i+1} = {round(result[i],round_num)}\n")
    print(f"Число итераций: {iters}")