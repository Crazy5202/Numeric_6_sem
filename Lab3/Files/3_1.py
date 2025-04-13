import os
import math

class NUMERIC_APPROX:
    def __init__(self):
        self.arg = 0.8
        self.a = []
        self.b = []

    def read_from_file(self, name):
        """Ввести информацию из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()

            self.arg = float(lines[0])
            
            line_stripped = lines[1].strip().split(' ')
            for elem in line_stripped:
                self.a.append(float(elem))

            line_stripped = lines[2].strip().split(' ')
            for elem in line_stripped:
                self.b.append(float(elem))

    def calc_eq(self, x):
        return math.e**x + x
    
    def find_index(self, array):
        index = 0
        for elem in array:
            if elem < self.arg:
                index += 1
            else:
                break
        return index

    def calc_lagr(self, order):
        sum = []
        for arr in [self.a, self.b]:
            if (order > len(arr)-1):
                print(f"\nПОРЯДОК СЛИШКОМ ВЫСОКИЙ! МАКСИМУМ {len(arr)-1}")
                return
            cur_sum = 0
            pos = self.find_index(arr)
            pos = min(pos, len(arr)-order-1)
            for i in range (pos, pos+order+1):
                cur_val = self.calc_eq(arr[i])
                for j in range (pos, pos+order+1):
                    if i==j:
                        continue
                    cur_val /= arr[i] - arr[j]
                    cur_val *= self.arg - arr[j]
                cur_sum += cur_val
            sum.append(cur_sum)
        print(f"\nЛагранж степени {order}:")
        for i in range (len(sum)):
            print(f"Погрешность для значений {i+1}: {abs(sum[i] - self.calc_eq(self.arg))}")

    def razd_razn_iter(self, array, left, right):
        if (left==right): return self.calc_eq(array[left])
        left_val = self.razd_razn_iter(array, left, right-1)
        right_val = self.razd_razn_iter(array, left+1, right)
        return (left_val - right_val) / (array[left] - array[right])
    
    def calc_newton(self, order):
        sum = []
        for arr in [self.a, self.b]:
            if (order > len(arr)-1):
                print(f"\nПОРЯДОК СЛИШКОМ ВЫСОКИЙ! МАКСИМУМ {len(arr)-1}")
                return
            cur_sum = 0
            pos = self.find_index(arr)
            pos = min(pos, len(arr)-order-1)
            for i in range (pos, pos+order+1):
                cur_val = self.razd_razn_iter(arr, pos, i)
                for j in range (pos, i):
                    cur_val *= (self.arg - arr[j])
                cur_sum += cur_val
            sum.append(cur_sum)
        print(f"\nНьютон степени {order}:")
        for i in range (len(sum)):
            print(f"Погрешность для значений {i+1}: {abs(sum[i] - self.calc_eq(self.arg))}")

if __name__ == "__main__":
    solver = NUMERIC_APPROX()

    solver.read_from_file("input_1.txt")
    solver.calc_lagr(2)
    solver.calc_lagr(3)
    solver.calc_newton(2)
    solver.calc_newton(3)