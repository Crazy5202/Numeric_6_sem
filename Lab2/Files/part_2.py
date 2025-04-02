import math
from matplotlib import pyplot as plt

class NONLINEAR_SLAU_SOLVER:
    def __init__(self):
        self.precision = 1e-8

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

    
    
    def draw_equation_1(self, x):
        return math.cos(x)/3
    
    def draw_equation_2(self, x):
        return math.log(3*x)
    
    def draw_equation(self):
        split = 1001
        left_edge = 0.1
        right_edge = 3
        
        x = [left_edge + (right_edge-left_edge)/split*i for i in range (split)]

        values = [self.draw_equation_1(x[i]) for i in range(len(x))]
        plt.plot(x, values)

        plt.vlines([left_edge], min(values), max(values), '0')

        values = [self.draw_equation_2(x[i]) for i in range(len(x))]
        plt.plot(x, values)
        
        plt.vlines([left_edge], min(values), max(values), '0')

        plt.hlines([0], left_edge, right_edge, '0')

        plt.title("Графики функции для выбора начального приближения.")
        plt.grid()
        plt.show()

    def calc_equation_1(self, x1, x2):
        return 3*x1 - math.cos(x2)
    
    def calc_equation_2(self, x1, x2):
        return 3*x2 - math.e**x1
    
    #
    #
    #

    def chosen_function_1(self, x1, x2):
        return math.cos(x2)/3

    def chosen_function_2(self, x1, x2):
        return math.e**x1/3
    
    def chosen_derivative_1(self, x1, x2):
        return -math.sin(x2)/3
    
    def chosen_derivative_2(self, x1, x2):
        return math.e**x1/3
    
    def check_conditions_iters(self, x1, x2, radius) -> tuple[bool, float]:
        max_coeff = 0
        for x in [x1-radius, x1, x1+radius]:
            for y in [x2-radius, x2, x2+radius]:
                matrix_norm = abs(self.chosen_derivative_1(x1,x2)) + abs(self.chosen_derivative_2(x1,x2))
                max_coeff = max(max_coeff, matrix_norm)
        if max_coeff > 1:
            print("Не выполняется ограниченность производной для данной функции и области.")
            return False, 0
        return True, max_coeff
    

    def simple_iters(self) -> tuple[float, int]:
        print("\n\nВыберите начальное приближение\n")
        ans = [0,0]
        counter = 0
        coeff = 0
        a = 0
        b = 0
        while(1):
            self.draw_equation()
        
            b = float(input("\nВведите первую координату: "))
            a = float(input("\nВведите вторую координату: "))
            r = float(input("\nВведите радиус области: "))
            fact, q = self.check_conditions_iters(a,b,r)
            if (fact):
                coeff = q/(1-q)
                ans = [a, b]
                #print(f"GOOD! {q}")
                break
        
        while (1):
            counter += 1
            new_ans = [self.chosen_function_1(ans[0], ans[1]), self.chosen_function_2(ans[0], ans[1])]
            if abs(new_ans[0] - a) > r or abs(new_ans[1] - b) > r:
                print("\n\nКорень не содержится в выбранной области!!!\n")
                return [], 0
            if (coeff * (abs(ans[0]-new_ans[0]) + abs(ans[1]-new_ans[1])) < self.precision):
                ans = new_ans
                break
            ans = new_ans
        return ans, counter
        
"""
    def check_conditions_newton(self, a, b) -> tuple[bool, float]:
        val = 0
        fact = True
        if self.calc_equation(a)*self.calc_equation(b)<0 and abs(self.calc_equation(a)*self.second_derivative(a))<self.first_derivative(a)**2 and (abs(self.calc_equation(b)*self.second_derivative(b))<self.first_derivative(b)**2):
            if (self.calc_equation(a)*self.second_derivative(a)>0):
                val = a
            elif (self.calc_equation(b)*self.second_derivative(b)>0):
                val = b
            else:
                fact = False
        else:
            fact = False
        if (fact == False):
            print("Границы выбраны неверно")
        return fact, val
        
    def newton(self) -> tuple[float,int]:
        """"Найти корень методом Ньютона.""""
        counter = 0
        ans = 0
        print("\n\nВыберите положительные границы (значения на границах должны различаться знаком)\n")
        while(1):
            self.draw_equation()
            a = float(input("\nВведите левую границу: "))
            b = float(input("\nВведите правую границу: "))
            if not(a>=0 and b>=a):
                print("Границы должны быть неотрицательными и вторая больше первой.")
                continue
            fact, ans = self.check_conditions_newton(a,b)
            if (fact):
                break
        while (1):
            counter += 1
            new_ans = ans - self.calc_equation(ans)/self.first_derivative(ans)
            if (abs(new_ans-ans)<self.precision):
                ans = new_ans
                break
            ans = new_ans
        return ans, counter 
"""

if __name__ == "__main__":
    solver = NONLINEAR_SLAU_SOLVER()

    round_num = solver.get_precision_num()
    
    #ans, iters = solver.newton()
    #print("\nМЕТОД НЬЮТОНА")
    #print(f"Корень: {ans}")

    #print(f"\nЧисло итераций: {iters}\n")
    
    ans, iters = solver.simple_iters()
    print("\nМЕТОД ПРОСТЫХ ИТЕРАЦИЙ")
    print(f"Корень: {ans}")

    print(f"Число итераций: {iters}")