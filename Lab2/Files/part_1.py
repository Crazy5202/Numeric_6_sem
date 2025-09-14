import math
from matplotlib import pyplot as plt

class NONLINEAR_SIGNLE_SOLVER:
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
    
    def second_derivative(self, x):
        return math.log(4)**2 * 4**x
    
    def first_derivative(self, x):
        return math.log(4) * 4**x - 5

    def calc_equation(self, x):
        return 4**x - 5*x - 2
    
    def draw_equation(self):
        split = 1001
        left_edge = -5
        right_edge = 5
        x = [left_edge + (right_edge-left_edge)/split*i for i in range (split)]
        values = [self.calc_equation(elem) for elem in x]
        plt.plot(x, [self.calc_equation(elem) for elem in x])
        plt.hlines([0], left_edge, right_edge, '0')
        plt.vlines([0], min(values), max(values), '0')
        plt.title("График функции для выбора начального приближения.")
        plt.grid()
        plt.show()

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
        """Найти корень методом Ньютона."""
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
    
             

    def chosen_derivative(self, x):
        return (5/(2*math.log(2)*(5*x+2)))

    def chosen_function(self, x):
        return math.log(5*x+2, 2)/2

    def check_conditions_iters(self, a, b) -> tuple[bool, float]:
        left_edge = a
        right_edge = b
        for elem in [left_edge, right_edge]:
            val = self.chosen_function(elem)
            if val<left_edge or val>right_edge:
                print("Функция не удовлетворяет условию 1 на отрезке")
                return False, 0
        if self.chosen_derivative(left_edge) >= 1 and self.chosen_derivative(right_edge)>=1:
            print("Функция не удовлетворяет условию 2 на отрезке")
            return False, 0
        max_val = min(self.chosen_derivative(left_edge), self.chosen_derivative(right_edge))
        return True, max_val

    def simple_iters(self) -> tuple[float, int]:
        print("\n\nВыберите положительные границы (значения на границах должны различаться знаком)\n")
        ans = 0
        coeff = 0
        counter = 0
        while(1):
            self.draw_equation()
            a = float(input("\nВведите левую границу: "))
            b = float(input("\nВведите правую границу: "))
            if not(a>=0 and b>=a):
                print("Границы должны быть неотрицательными и вторая больше первой.")
                continue
            fact, q = self.check_conditions_iters(a,b)
            if (fact):
                ans = b
                coeff = q/(1-q)
                break
        while (1):
            counter += 1
            new_ans = self.chosen_function(ans)
            if (coeff*abs(new_ans-ans)<self.precision):
                ans = new_ans
                break
            ans = new_ans
        return ans, counter

if __name__ == "__main__":
    solver = NONLINEAR_SIGNLE_SOLVER()

    round_num = solver.get_precision_num()
    
    ans, iters = solver.newton()
    print("\nМЕТОД НЬЮТОНА")
    print(f"Корень: {ans}")
    print(f"\nЧисло итераций: {iters}\n")
    
    ans, iters = solver.simple_iters()
    print("\nМЕТОД ПРОСТЫХ ИТЕРАЦИЙ")
    print(f"Корень: {ans}")
    print(f"\nЧисло итераций: {iters}\n")
    