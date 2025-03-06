import os

class TriDiag:
    def __init__(self):
        self.a = []
        self.b = []
        self.c = []
        self.d = []
        self.n = 0

    def read_from_file(self, name):
        """Ввести матрицу из файла"""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()
            self.n = len(lines)

            line_stripped = lines[0].strip().split(' ')
            self.a.append(0)
            self.b.append(float(line_stripped[0]))
            self.c.append(float(line_stripped[1]))
            self.d.append(float(line_stripped[2]))

            for i in range(1,self.n-1):
                line_stripped = lines[i].strip().split(' ')
                self.a.append(float(line_stripped[0]))
                self.b.append(float(line_stripped[1]))
                self.c.append(float(line_stripped[2]))
                self.d.append(float(line_stripped[3]))

            line_stripped = lines[self.n-1].strip().split(' ')
            self.a.append(float(line_stripped[0]))
            self.b.append(float(line_stripped[1]))
            self.c.append(0)
            self.d.append(float(line_stripped[2]))

    def check_conditions(self):
        """Проверяет корректность и устойчивость."""
        for i in range (self.n):
            if not (abs(self.b[i]) > 0 and abs(self.b[i]) >= abs(self.a[i]) + abs(self.c[i])):
                return False
        return True
    
    def check_solution(self, x):
        calc_d = []
        calc_d.append(self.b[0]*x[0]+self.c[0]*x[1])
        for i in range(1,self.n-1):
            calc_d.append(self.a[i]*x[i-1]+self.b[i]*x[i]+self.c[i]*x[i+1])
        calc_d.append(self.a[self.n-1]*x[self.n-2]+self.b[self.n-1]*x[self.n-1])
        for i in range (self.n):
            diff = calc_d[i] - self.d[i]
            if (diff > 1e-6):
                return False
        return True 

    def solve(self):
        """Решает методом прогонки."""
        if self.check_conditions() == False:
            print("Не выполняется проверка!")
            return []
        A = []
        B = []
        x = [0 for _ in range(self.n)]
        A.append(-self.c[0]/self.b[0])
        B.append(self.d[0]/self.b[0])
        for i in range (1,self.n):
            A.append(-self.c[i]/(self.b[i]+self.a[i]*A[i-1]))
            B.append((self.d[i]-self.a[i]*B[i-1])/(self.b[i]+self.a[i]*A[i-1]))
        x[self.n-1] = B[self.n-1]
        for i in range (self.n-2, -1, -1):
            x[i] = A[i]*x[i+1]+B[i]
        if (self.check_solution(x)): 
            print("\nРешение верно!\n")
            return x
        else:
            print("\nРешение неверно!\n")
            return []

    def print_matrix(self):
        """Печатает матрицу."""
        for i in range (self.n):
            print(self.a[i], self.b[i], self.c[i], "=", self.d[i], "\n")

if __name__ == "__main__":
    matrix = TriDiag()
    matrix.read_from_file("matrix_2.txt")
    result = matrix.solve()
    for i in range (len(result)):
        print(f"x_{i} = {round(result[i],3)}\n")