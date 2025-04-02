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
    
    def eq1_d1(self, x1, x2):
        return 3
    
    def eq1_d2(self, x1, x2):
        return math.sin(x2)
    
    def eq2_d1(self, x1, x2):
        return -math.e**x1
    
    def eq2_d2(self, x1, x2):
        return 3
    
    def jacobian_norm(self, x1, x2):
        return self.eq1_d1(x1,x2)*self.eq2_d2(x1,x2) - self.eq1_d2(x1,x2)*self.eq2_d1(x1,x2)
        
    def calc_d(self, x1, x2):
        matrix = [[self.eq1_d1(x1,x2), self.eq1_d2(x1,x2), self.calc_equation_1(x1,x2)], 
                  [self.eq2_d1(x1,x2), self.eq2_d2(x1,x2), self.calc_equation_2(x1,x2)]]
        coeff = matrix[1][0]/matrix[0][0]
        for i in range (len(matrix[0])):
            matrix[1][i] -= matrix[0][i]
        coeff = matrix[0][1]/matrix[1][1]
        for i in range (len(matrix[0])):
            matrix[0][i] -= matrix[1][i]
        return [matrix[0][2]/matrix[0][0], matrix[1][2]/matrix[1][1]]


    def newton(self) -> tuple[float,int]:
        """Найти корень методом Ньютона."""
        counter = 0
        a = 0
        b = 0
        print("\n\nВыберите начальное приближение\n")
        self.draw_equation()
        b = float(input("\nВведите первую координату: "))
        a = float(input("\nВведите вторую координату: "))
        
        ans = [a, b]

        while (1):
            counter += 1
            if self.jacobian_norm(ans[0], ans[1]) < 1e-6:
                return [], 0
            d = self.calc_d(ans[0], ans[1])
            new_ans = [ans[0]-d[0], ans[1]-d[1]] 
            if (abs(new_ans[0]-ans[0])+abs(new_ans[1]-ans[1])<self.precision):
                ans = new_ans
                break
            ans = new_ans
        return ans, counter 
    



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
        """Найти корень методом простых итераций."""
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




if __name__ == "__main__":
    solver = NONLINEAR_SLAU_SOLVER()

    round_num = solver.get_precision_num()
    
    ans, iters = solver.newton()
    print("\nМЕТОД НЬЮТОНА")
    print(f"Корень: {round(ans[0], round_num)}, {round(ans[1], round_num)}")
    print(f"\nЧисло итераций: {iters}\n")
    
    
    ans, iters = solver.simple_iters()
    print("\nМЕТОД ПРОСТЫХ ИТЕРАЦИЙ")
    print(f"Корень: {round(ans[0], round_num)}, {round(ans[1], round_num)}")
    print(f"Число итераций: {iters}")
    