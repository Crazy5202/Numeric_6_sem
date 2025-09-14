import os
import math

class NUMERIC_INTEGR:
    def __init__(self):
        self.x0 = 0
        self.x1 = 0
        self.h1 = 0
        self.h2 = 0

    def read_from_file(self, name):
        """Ввести информацию из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()

            line_stripped = lines[0].strip().split(' ')
            self.x0 = float(line_stripped[0])
            self.x1 = float(line_stripped[1])

            line_stripped = lines[1].strip().split(' ')
            self.h1 = float(line_stripped[0])
            self.h2 = float(line_stripped[1])
        print("DONE!")

    def calc_eq(self, x):
        return 1 / (256 - x**4)
    
    def runge(self, sum1, sum2):
        actual_value = 1/128*math.log(3) + math.atan(1/2)/64
        refined = sum1 + (sum1 - sum2)/((self.h2/self.h1)**2 - 1)
        print(f"С шагом {self.h1}: {sum1}, с шагом {self.h2}: {sum2}")
        print(f"Уточнённый по Рунге: {refined}, погрешность: {abs(actual_value-refined)}")

    def calc_rect(self):
        sum = []
        for h in [self.h1, self.h2]:
            cur_sum = 0
            x = self.x0
            while (x < self.x1):
                cur_sum += self.calc_eq(x + h/2)
                x += h
            sum.append(h * cur_sum)
        print("\n\nМЕТОД ПРЯМОУГОЛЬНИКОВ")
        self.runge(sum[0], sum[1])

    def calc_trap(self):
        sum = []
        for h in [self.h1, self.h2]:
            cur_sum = 0
            x = self.x0
            while (x < self.x1):
                cur_sum += self.calc_eq(x) + self.calc_eq(x+h)
                x += h
            sum.append(h * cur_sum / 2)
        print("\n\nМЕТОД ТРАПЕЦИЙ")
        self.runge(sum[0], sum[1])

    def calc_simp(self):
        sum = []
        for h in [self.h1, self.h2]:
            cur_sum = 0
            x = self.x0
            while (x < self.x1):
                cur_sum +=  self.calc_eq(x) + 4*self.calc_eq(x+h/2) + self.calc_eq(x+h)
                x += h
            sum.append(h/2 * cur_sum / 3)
        print("\n\nМЕТОД СИМПСОНА")
        self.runge(sum[0], sum[1])

if __name__ == "__main__":
    solver = NUMERIC_INTEGR()

    solver.read_from_file("input_5.txt")
    solver.calc_rect()
    solver.calc_trap()
    solver.calc_simp()