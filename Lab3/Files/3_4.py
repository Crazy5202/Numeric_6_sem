import os

class NUMERIC_DIFF:
    def __init__(self):
        self.arg = 0
        self.x = []
        self.y = []

    def read_from_file(self, name):
        """Ввести информацию из файла."""
        PATH = os.path.split(os.path.realpath(__file__))[0] + "/" + name
        with open(PATH, "r") as f:
            lines = f.readlines()
            self.arg = float(lines[0])
            for i in range(1, len(lines)):
                line_stripped = lines[i].strip().split(' ')
                self.x.append(float(line_stripped[0]))
                self.y.append(float(line_stripped[1]))

    def find_index(self):
        index = 0
        for elem in self.x:
            if elem < self.arg:
                index += 1
            else:
                break
        return index
    
    def left_diff(self, ind):
        return (self.y[ind]-self.y[ind-1])/(self.x[ind]-self.x[ind-1])

    def right_diff(self, ind):
        return (self.y[ind+1]-self.y[ind])/(self.x[ind+1]-self.x[ind])

    def first_diff(self) -> float:
        """Посчитать первую производную."""
        index = self.find_index()

        return (self.left_diff(index) + self.right_diff(index))/2

    def second_diff(self) -> float:
        """Посчитать вторую производную."""

        index = self.find_index()

        return 2*((self.right_diff(index) - self.left_diff(index))) / (self.x[index+1] - self.x[index-1])
                

if __name__ == "__main__":
    solver = NUMERIC_DIFF()

    solver.read_from_file("input_4.txt")
    print(f"Первая производная: {solver.first_diff()}")
    print(f"Вторая производная: {solver.second_diff()}")