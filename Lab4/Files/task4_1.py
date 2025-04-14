from typing import Callable
import math
from multiprocessing import Process


def linspace(a: float, b: float, points: float) -> list[float]:
    return [a + (b - a) / (points - 1) * i for i in range(points)]


class ODE2_Euler:
    def __init__(self, x0: float, xn: float, h: float, y0: float, z0: float, f: Callable, g: Callable):
        X = linspace(x0, xn, int((xn - x0) / h + 1))
        Y = [y0 for _ in range(len(X))]
        Z = [z0 for _ in range(len(X))]

        for i in range(1, len(X)):
            Z[i] = Z[i - 1] + h * g(X[i - 1], Y[i - 1], Z[i - 1])
            Y[i] = Y[i - 1] + h * f(X[i - 1], Y[i - 1], Z[i - 1])

        self.X = X
        self.Y = Y
        self.Z = Z

    def __call__(self, x: float) -> float:
        if x < self.X[0] or x > self.X[-1]:
            raise Exception('X is out of range')

        for i in range(1, len(self.X)):
            if self.X[i - 1] <= x <= self.X[i]:
                return self.Y[i - 1] + (self.Y[i] - self.Y[i - 1]) * (x - self.X[i - 1]) / (self.X[i] - self.X[i - 1])

    def MAE(self, true: Callable) -> float:
        return sum([math.fabs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])


class ODE2_Runge_Kutta:
    def __init__(self, x0: float, xn: float, h: float, y0: float, z0: float, f: Callable, g: Callable):
        X = linspace(x0, xn, int((xn - x0) / h + 1))
        Y = [y0 for _ in range(len(X))]
        Z = [z0 for _ in range(len(X))]

        for i in range(len(X) - 1):
            K1 = h * f(X[i], Y[i], Z[i])
            L1 = h * g(X[i], Y[i], Z[i])

            K2 = h * f(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K1, Z[i] + 1 / 2 * L1)
            L2 = h * g(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K1, Z[i] + 1 / 2 * L1)

            K3 = h * f(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K2, Z[i] + 1 / 2 * L2)
            L3 = h * g(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K2, Z[i] + 1 / 2 * L2)

            K4 = h * f(X[i] + h, Y[i] + K3, Z[i] + L3)
            L4 = h * g(X[i] + h, Y[i] + K3, Z[i] + L3)

            Y[i + 1] = Y[i] + 1 / 6 * (K1 + 2 * K2 + 2 * K3 + K4)
            Z[i + 1] = Z[i] + 1 / 6 * (L1 + 2 * L2 + 2 * L3 + L4)

        self.X = X
        self.Y = Y
        self.Z = Z

    def __call__(self, x: float) -> float:
        if x < self.X[0] or x > self.X[-1]:
            raise Exception('X is out of range')

        for i in range(1, len(self.X)):
            if self.X[i - 1] <= x <= self.X[i]:
                return self.Y[i - 1] + (self.Y[i] - self.Y[i - 1]) * (x - self.X[i - 1]) / (self.X[i] - self.X[i - 1])

    def MAE(self, true: Callable) -> float:
        return sum([math.fabs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])


class ODE2_Adams:
    def __init__(self, x0: float, xn: float, h: float, y0: float, z0: float, f: Callable, g: Callable):
        X = linspace(x0, xn, int((xn - x0) / h + 1))
        Y = [y0 for _ in range(len(X))]
        Z = [z0 for _ in range(len(X))]
        F = [f(x0, y0, z0) for _ in range(len(X))]
        G = [g(x0, y0, z0) for _ in range(len(X))]

        for i in range(3):
            K1 = h * f(X[i], Y[i], Z[i])
            L1 = h * g(X[i], Y[i], Z[i])

            K2 = h * f(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K1, Z[i] + 1 / 2 * L1)
            L2 = h * g(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K1, Z[i] + 1 / 2 * L1)

            K3 = h * f(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K2, Z[i] + 1 / 2 * L2)
            L3 = h * g(X[i] + 1 / 2 * h, Y[i] + 1 /
                       2 * K2, Z[i] + 1 / 2 * L2)

            K4 = h * f(X[i] + h, Y[i] + K3, Z[i] + L3)
            L4 = h * g(X[i] + h, Y[i] + K3, Z[i] + L3)

            Y[i + 1] = Y[i] + 1 / 6 * (K1 + 2 * K2 + 2 * K3 + K4)
            Z[i + 1] = Z[i] + 1 / 6 * (L1 + 2 * L2 + 2 * L3 + L4)
            F[i + 1] = f(X[i + 1], Y[i + 1], Z[i + 1])
            G[i + 1] = g(X[i + 1], Y[i + 1], Z[i + 1])

        for i in range(3, len(X) - 1):
            Y[i + 1] = Y[i] + h / 24 * \
                (55 * F[i] - 59 * F[i - 1] + 37 * F[i - 2] - 9 * F[i - 3])

            Z[i + 1] = Z[i] + h / 24 * \
                (55 * G[i] - 59 * G[i - 1] + 37 * G[i - 2] - 9 * G[i - 3])

            F[i + 1] = f(X[i + 1], Y[i + 1], Z[i + 1])
            G[i + 1] = g(X[i + 1], Y[i + 1], Z[i + 1])

        self.X = X
        self.Y = Y
        self.Z = Z

    def __call__(self, x: float) -> float:
        if x < self.X[0] or x > self.X[-1]:
            raise Exception('X is out of range')

        for i in range(1, len(self.X)):
            if self.X[i - 1] <= x <= self.X[i]:
                return self.Y[i - 1] + (self.Y[i] - self.Y[i - 1]) * (x - self.X[i - 1]) / (self.X[i] - self.X[i - 1])

    def MAE(self, true: Callable) -> float:
        return sum([math.fabs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])


class RungeRomberg:
    def __init__(self, X: list[float], Y_h: list[float], Y_h2: list[float], p: float):
        self.X = X
        self.Y = [y_h2 + (y_h2 - y_h) / (2**p - 1)
                  for y_h, y_h2 in zip(Y_h, Y_h2[::2])]

    def __call__(self, x: float) -> float:
        if x < self.X[0] or x > self.X[-1]:
            raise Exception('X is out of range')

        for i in range(1, len(self.X)):
            if self.X[i - 1] <= x <= self.X[i]:
                return self.Y[i - 1] + (self.Y[i] - self.Y[i - 1]) * (x - self.X[i - 1]) / (self.X[i] - self.X[i - 1])

    def MAE(self, true: Callable) -> float:
        return sum([math.fabs(self.Y[i] - true(self.X[i])) for i in range(len(self.X))])


def f(x: float, y: float, z: float) -> float:
    return z


def g(x: float, y: float, z: float) -> float:
    return 12 * y / x**2


def y(x: float) -> float:
    return x**4 + 1 / x**3


def main():
    x0 = 1
    xn = 2
    h = 0.1
    y0 = 2
    z0 = 1

    euler1 = ODE2_Euler(x0, xn, h, y0, z0, f, g)
    euler2 = ODE2_Euler(x0, xn, h / 2, y0, z0, f, g)
    rr_euler = RungeRomberg(euler1.X, euler1.Y, euler2.Y, 2)

    runge1 = ODE2_Runge_Kutta(x0, xn, h, y0, z0, f, g)
    runge2 = ODE2_Runge_Kutta(x0, xn, h / 2, y0, z0, f, g)
    rr_runge = RungeRomberg(runge1.X, runge1.Y, runge2.Y, 4)

    adams1 = ODE2_Adams(x0, xn, h, y0, z0, f, g)
    adams2 = ODE2_Adams(x0, xn, h / 2, y0, z0, f, g)
    rr_adams = RungeRomberg(adams1.X, adams1.Y, adams2.Y, 4)

    print('Метод Эйлера: ')
    print(f'Абсолютная погрешность: {euler1.MAE(y)}')
    print(
        f'Погрешность при применении метода Рунге-Ромберга: {rr_euler.MAE(y)}', end='\n\n')

    print('Метод Рунге-Кутты: ')
    print(f'Абсолютная погрешность: {runge1.MAE(y)}')
    print(
        f'Погрешность при применении метода Рунге-Ромберга: {rr_runge.MAE(y)}', end='\n\n')

    print('Метод Адамса: ')
    print(f'Абсолютная погрешность: {adams1.MAE(y)}')
    print(
        f'Погрешность при применении метода Рунге-Ромберга: {rr_adams.MAE(y)}')


if __name__ == "__main__":
    main()
