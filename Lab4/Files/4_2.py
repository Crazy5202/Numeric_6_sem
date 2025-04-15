import os
import math
from matplotlib import pyplot as plt

class NUMERIC_KOSHI_23:
    def __init__(self):
        """Параметры: x0, x1, первое условие в x0, второе условие в x0."""
        self.x0 = 0
        self.x1 = 1
        self.corner1 = 0
        self.precision = 1e-3
    
    def calc_eq(self, x, y, z):
        return (y-(x-3)*z)/(x**2-1)
    
    def calc_algor(self, z, y):
        val = z + y + 0.75
        return val

    def calc_solution(self, x):
        return x - 3 + 1/(x+1)
    
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
        """Уточнить значение по Рунге-Ромбергу."""
        refined = [sum1[i] + (sum1[i] - sum2[i*2])/((1/2)**2 - 1) for i in range(len(sum1))]
        return refined

    def calc_kutta(self, h, corner0) -> tuple[list, list]:
        """Рассчитать значения по Рунге-Кутте."""
        x_prev = self.x0
        y_prev = corner0
        z_prev = self.corner1

        y_array = [y_prev]
        z_array = [z_prev]
        
        while (x_prev <= self.x1-h/10):
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

            z_prev = z_new
            y_prev = y_new

            y_array.append(y_prev)
            z_array.append(z_prev)
            
            x_prev += h

        return z_array, y_array
    
    def calc_eiler(self, h, corner0):
        """Рассчитать значения по Эйлеру."""
        
        x_prev = self.x0
        y_prev = corner0
        z_prev = self.corner1

        y_array = [y_prev]
        z_array = [z_prev]

        while (x_prev < self.x1 - h/10):
            y_new = y_prev + h * z_prev
            z_new = z_prev + h * self.calc_eq(x_prev, y_prev, z_prev)

            y_array.append(y_new)
            z_array.append(z_new)

            z_prev = z_new
            y_prev = y_new
            
            x_prev += h

        return z_array, y_array
    
    def calc_shooting(self, h) -> list:
        """Рассчитать значения методом стрельбы."""
        method = self.calc_eiler

        arg0 = 0.3
        arg1 = 0.5
        
        res0 = method(h, arg0)
        res1 = method(h, arg1)

        algor0 = self.calc_algor(res0[0][-1], res0[1][-1])
        algor1 = self.calc_algor(res1[0][-1], res1[1][-1])
        
        while abs(algor1) > self.precision:
            new_arg = arg1 - (arg1-arg0)/(algor1-algor0)*algor1
            new_res = method(h, new_arg)
            new_algor = self.calc_algor(new_res[0][-1], new_res[1][-1])

            arg0 = arg1
            arg1 = new_arg
            
            algor0 = algor1
            algor1 = new_algor
        return new_res[1]
    
    def wrapper(self, h, method):
        """Обернуть метод для выполнения общего набора действий."""
        sum1 = method(h)
        sum2 = method(h/2)
        refined = self.runge(sum1, sum2)
        x_array = [self.x0 + h*i for i in range(int((self.x1-self.x0)/h)+1)]
        
        print(f"Абсолютная ошибка для сетки {h}: {self.calc_error(x_array, sum1)}")
        print(f"Абсолютная ошибка с уточнением по Рунге-Ромбергу: {self.calc_error(x_array, refined)}\n\n")

        self.plot(x_array, sum1, refined)

    def wrapper_wrapper(self, h):
        """Вызвать все методы."""
        print("\n\nМЕТОД СТРЕЛЬБЫ\n")
        self.wrapper(h, self.calc_shooting)

if __name__ == "__main__":
    solver = NUMERIC_KOSHI_23()

    solver.wrapper_wrapper(0.01)