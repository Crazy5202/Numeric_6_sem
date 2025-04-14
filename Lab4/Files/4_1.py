import os
import math
from matplotlib import pyplot as plt

class NUMERIC_KOSHI:
    def __init__(self, x0, x1, corner0, corner1):
        """Параметры: x0, x1, первое условие в x0, второе условие в x0."""
        self.x0 = x0
        self.x1 = x1
        self.corner0 = corner0
        self.corner1 = corner1
    
    def calc_eq(self, x, y, z):
        return (1/x+1)*z - y/x

    def calc_solution(self, x):
        return x + 1 + math.e**x
    
    def plot(self, x_array, integrated, refined):
        plt.plot(x_array, [solver.calc_solution(elem) for elem in x_array], label = 'Правильное решение')
        plt.plot(x_array, integrated, label = 'Численное решение')
        plt.plot(x_array, refined, label = 'Уточнённое численное решение')
        plt.legend()
        plt.grid()
        plt.show()
    
    def calc_error(self, x_array, results):
        sum = 0
        
        for i in range (len(results)):
            sum += abs(results[i] - self.calc_solution(x_array[i]))

        return sum
    
    def runge(self, sum1, sum2):
        """Уточнить значение по Рунге."""
        refined = [sum1[i] + (sum1[i] - sum2[i*2+1])/((1/2)**2 - 1) for i in range(len(sum1))]
        return refined

    def calc_eiler(self, h) -> list:
        """Рассчитать значения по Эйлеру."""
        result = []
        x_prev = self.x0
        y_prev = self.corner0
        z_prev = self.corner1
        while (x_prev < self.x1):
            y_new = y_prev + h * z_prev
            z_new = z_prev + h * self.calc_eq(x_prev, y_prev, z_prev)

            result.append(y_new)

            z_prev = z_new
            y_prev = y_new
            
            x_prev += h

        return result

    def eiler_wrapper(self, h):
        sum1 = self.calc_eiler(h)
        sum2 = self.calc_eiler(h/2)
        refined = self.runge(sum1, sum2)
        x_array = [self.x0 + h*(i+1) for i in range(int((self.x1-self.x0)/h))]
        
        print("\n\nМЕТОД ЭЙЛЕРА\n")
        print(f"Абсолютная ошибка для сетки {h}: {self.calc_error(x_array, sum1)}")
        print(f"Абсолютная ошибка с уточнением по Рунге-Ромбергу: {self.calc_error(x_array, refined)}\n\n")

        self.plot(x_array, sum1, refined)
        return refined
    
    def calc_kutta(self, h):
        
        result = []
        x_prev = self.x0
        y_prev = self.corner0
        z_prev = self.corner1
        while (x_prev < self.x1):
            y1 = h * z_prev
            z1 = h * self.calc_eq(x_prev, y_prev, z_prev)

            y2 = h * (z_prev + 1/2*z1)
            z2 = h * self.calc_eq(x_prev + 1/2*h, y_prev + 1/2*y1, z_prev + 1/2*z1)

            y3 = h * (z_prev + 1/2*z2)
            z3 = h * self.calc_eq(x_prev + 1/2*h, y_prev + 1/2*y2, z_prev + 1/2*z2)

            y4 = h * (z_prev + z3)
            z4 = h * self.calc_eq(x_prev + h, y_prev + y3, z_prev + z3)

            y_new = y_prev + 1/6*(y1 + 2*y2 + 2*y3 + y4)
            z_new = z_prev + 1/6*(z1 + 2*z2 + 2*z3 + z4)

            result.append(y_new)

            z_prev = z_new
            y_prev = y_new
            
            x_prev += h

        return result
    
    def kutta_wrapper(self, h):
        sum1 = self.calc_kutta(h)
        sum2 = self.calc_kutta(h/2)
        refined = self.runge(sum1, sum2)
        x_array = [self.x0 + h*(i+1) for i in range(int((self.x1-self.x0)/h))]
        
        print("\n\nМЕТОД РУНГЕ-КУТТЫ\n")
        print(f"Абсолютная ошибка для сетки {h}: {self.calc_error(x_array, sum1)}")
        print(f"Абсолютная ошибка с уточнением по Рунге-Ромбергу: {self.calc_error(x_array, refined)}\n\n")

        self.plot(x_array, sum1, refined)

if __name__ == "__main__":
    solver = NUMERIC_KOSHI(1, 2, 2 + math.e, 1 + math.e)

    h = 0.1

    #solver.eiler_wrapper(h)
    #solver.kutta_wrapper(h)