import os
import math
from matplotlib import pyplot as plt

from tridiag import TRIDIAG_SOLVER


class NUMERIC_APPROX_SPLINE:
    def __init__(self):
        self.arg = 0.0
        self.x = []
        self.f = []
        self.coeffs = []

    def read_from_file(self, name):
        """Ввести информацию из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()

            self.arg = float(lines[0])
            
            line_stripped = lines[1].strip().split(' ')
            for elem in line_stripped:
                self.x.append(float(elem))

            line_stripped = lines[2].strip().split(' ')
            for elem in line_stripped:
                self.f.append(float(elem))
    
    def find_index(self, array = None, argument = None):
        """Найти подходящий индекс для интерполяции."""
        if array is None:
            array = self.x
        if argument is None:
            argument = self.arg
        index = 0
        for elem in array:
            if elem < argument:
                index += 1
            else:
                break
        if index != 0:
            index -= 1
        return index

    def calc_h(self, ind):
        return self.x[ind] - self.x[ind-1]

    def prepare_spline_coeffs(self):
        solver = TRIDIAG_SOLVER()
        name = "/temp"
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/temp"

        n = len(self.x)-1

        with open(PATH, "w") as f:
            f.write(f"{ 2*(self.calc_h(1) + self.calc_h(2)) }")
            f.write(f" { self.calc_h(2) }")
            f.write(f" { 3*( (self.f[2]-self.f[1])/self.calc_h(2) - (self.f[1]-self.f[0])/self.calc_h(1) ) }\n")

            for i in range(3, n):
                f.write(f"{ self.calc_h(i-1) }")
                f.write(f" { 2*(self.calc_h(i-1) + self.calc_h(i)) }")
                f.write(f" { self.calc_h(i) }")
                f.write(f" { 3*( (self.f[i]-self.f[i-1])/self.calc_h(i) - (self.f[i-1]-self.f[i-2])/self.calc_h(i-1) ) }\n")

            f.write(f"{ self.calc_h(n-1) }")
            f.write(f" { 2*(self.calc_h(n-1) + self.calc_h(n)) }")
            f.write(f" { 3*( (self.f[n]-self.f[n-1])/self.calc_h(n) - (self.f[n-1]-self.f[n-2])/self.calc_h(n-1) ) }")
        solver.read_from_file("temp")
        c = [0] + solver.solve()

        os.remove(PATH)
        
        a = [self.f[i-1] for i in range(1,n+1)]
        for i in range(1, n):
            b = (self.f[i]-self.f[i-1])/self.calc_h(i) - 1/3*self.calc_h(i)*(c[i]+2*c[i-1])
            d = (c[i]-c[i-1])/(3*self.calc_h(i))
            self.coeffs.append([a[i-1], b, c[i-1], d])
        b = (self.f[n]-self.f[n-1])/self.calc_h(n) - 2/3*self.calc_h(n)*c[n-1]
        d = -c[n-1]/(3*self.calc_h(n))
        self.coeffs.append([a[n-1], b, c[-1],d])
    
    def calc_spline(self, x):
        ind = self.find_index(argument=x)
        x_prev = self.x[ind]
        coeffs = self.coeffs[ind]
        sum = 0
        for i in range(len(coeffs)):
            sum += coeffs[i]*(x - x_prev)**(i)
        return sum

if __name__ == "__main__":
    solver = NUMERIC_APPROX_SPLINE()

    solver.read_from_file("input_2.txt")
    solver.prepare_spline_coeffs()

    print(f"\nРешение для {solver.arg}: {solver.calc_spline(solver.arg)}\n\n")
    
    split = 100
    left_edge = solver.x[0]
    right_edge = solver.x[-1]
        
    x = [left_edge + (right_edge-left_edge)/split*i for i in range (split+1)]

    func_values = [solver.calc_spline(elem) for elem in x]
    plt.plot(x, func_values)
    plt.scatter(solver.x, solver.f)
    plt.title("Кубический сплайн")
    plt.grid()
    plt.show()