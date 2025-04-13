import os
import math
from matplotlib import pyplot as plt

from LU import LU_SOLVER

class NUMERIC_APPROX_MNC:
    def __init__(self):
        self.x = []
        self.y = []
        self.coeffs = []

    def read_from_file(self, name):
        """Ввести информацию из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()
            
            line_stripped = lines[0].strip().split(' ')
            for elem in line_stripped:
                self.x.append(float(elem))

            line_stripped = lines[1].strip().split(' ')
            for elem in line_stripped:
                self.y.append(float(elem))

    def prepare_coeffs(self, order):
        self.coeffs = []

        name = "temp"

        solver = LU_SOLVER()
        PATH = os.path.split(os.path.realpath(__file__))[0] + f"/{name}"

        with open(PATH, "w") as f:
            for k in range(order+1):
                for i in range(order+1):
                    coeff = 0
                    for j in range(len(self.x)):
                        coeff += self.x[j]**(k+i)
                    f.write(f"{coeff} ")
                coeff = 0
                for j in range (len(self.x)):
                    coeff += (self.x[j]**k) * self.y[j]
                f.write(f"{coeff}\n")

        solver.read_from_file(name)
        self.coeffs = solver.solve()

        os.remove(PATH)

        error = 0

        for i in range (len(self.x)):
            error += (self.calc_approx(self.x[i]) - self.y[i])**2

        print(f"\nДля порядка {order} квадратичное отклонение: {error}\n")
    
    def calc_approx(self, x):
        sum = 0
        for i in range(len(self.coeffs)):
            sum += self.coeffs[i]*(x**i)
        return sum

if __name__ == "__main__":
    solver = NUMERIC_APPROX_MNC()

    solver.read_from_file("input_3.txt")
    
    split = 100
    left_edge = solver.x[0]
    right_edge = solver.x[-1]
    
    x = [left_edge + (right_edge-left_edge)/split*i for i in range (split+1)]

    solver.prepare_coeffs(1)

    func_values = [solver.calc_approx(elem) for elem in x]
    plt.plot(x, func_values, label = '1 степень')

    solver.prepare_coeffs(2)

    func_values = [solver.calc_approx(elem) for elem in x]
    plt.plot(x, func_values, label = '2 степень')

    plt.scatter(solver.x, solver.y, color = 'black')
    plt.title("Метод наименьших квадратов")
    plt.grid()
    plt.legend()
    plt.show()
    