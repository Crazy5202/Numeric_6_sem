from task4_1 import ODE2_Runge_Kutta, RungeRomberg
import math
from typing import Callable
from multiprocessing import Process


def f(x: float, y: float, z: float) -> float:
    return z


def g(x: float, y: float, z: float) -> float:
    return (y-(x-3)*z)/(x**2-1)


def Phi(solution: ODE2_Runge_Kutta) -> float:
    return solution.Y[-1] + solution.Z[-1] + 0.75


def p(x: float) -> float:
    return 2 * (x + 1) / (x * (2 * x + 1))


def q(x: float) -> float:
    return - 2 / (x * (2 * x + 1))


def y(x: float) -> float:
    return x - 3 + 1/(x+1)


def shooting_method(x0: float, xn: float, h: float, z0: float, eps=1e-5) -> ODE2_Runge_Kutta:
    s0 = 1.0
    s1 = 0.8

    solution0 = ODE2_Runge_Kutta(x0, xn, h, s0, z0, f, g)
    phi0 = Phi(solution0)
    if math.fabs(phi0) < eps:
        return solution0

    solution1 = ODE2_Runge_Kutta(x0, xn, h, s1, z0, f, g)
    phi1 = Phi(solution1)

    while math.fabs(phi1) >= eps:
        sk = s1 - (s1 - s0) / (phi1 - phi0) * phi1
        s0, s1 = s1, sk

        solution1 = ODE2_Runge_Kutta(x0, xn, h, s1, z0, f, g)
        phi0 = phi1
        phi1 = Phi(solution1)

    return solution1



def main():
    x0 = 1
    xn = 2
    h = 0.1
    z0 = 0

    shooting1 = shooting_method(x0, xn, h, z0)
    shooting2 = shooting_method(x0, xn, h / 2, z0)
    rr_shooting = RungeRomberg(shooting1.X, shooting1.Y, shooting2.Y, 2)


    print('Метод стрельбы:')
    print(f'Абсолютная погрешность: {shooting1.MAE(y)}')
    print(
        f'Погрешность при использовании метода Рунге-Ромберга: {rr_shooting.MAE(y)}', end='\n\n')


if __name__ == "__main__":
    main()
