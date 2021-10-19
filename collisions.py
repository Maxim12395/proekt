import numpy as np

def collision(x1, y1, vx1, vy1, x2, y2, vx2, vy2, radius_1, radius_2, mass1, mass2, K):

    """Аргументы функции:
    x1,y1,vx1,vy1 - координаты и компоненты скорости 1-ой частицы
    x2,y2,vx2,vy2 - ... 2-ой частицы
    radius,mass1,mass2 - радиус частиц и их массы (массы разные можно задавать,
    радиус для простоты взят одинаковый)
    K - коэффициент восстановления (K=1 для абсолютного упругого удара, K=0
    для абсолютно неупругого удара, 0<K<1 для реального удара)
    Функция возвращает компоненты скоростей частиц, рассчитанные по формулам для
    реального удара, если стокновение произошло. Если удара нет, то возвращаются
    те же значения скоростей, что и заданные в качестве аргументов.
    """
    r12 = np.sqrt((x1-x2)**2 + (y1-y2)**2) #расчет расстояния между центрами частиц
    v1 = np.sqrt(vx1**2 + vy1**2) #расчет модулей скоростей частиц
    v2 = np.sqrt(vx2**2 + vy2**2)


    # проверка условия на столкновение: расстояние должно быть меньше 2-х радиусов
    if r12 <= radius_1 + radius_2:
        '''вычисление углов движения частиц theta1(2), т.е. углов между
        направлением скорости частицы и положительным направлением оси X.
        Если частица  покоится, то угол считается равным нулю. Т.к. функция
        arccos имеет область значений от 0 до Pi, то в случае отрицательных
        y-компонент скорости для вычисления угла theta1(2) надо из 2*Pi
        вычесть значение arccos(vx/v)
        '''
        if v1 != 0:
            theta1 = np.arccos(vx1 / v1)
        else:
            theta1 = 0
        if v2 != 0:
            theta2 = np.arccos(vx2 / v2)
        else:
            theta2 = 0
        if vy1 < 0:
            theta1 = - theta1 + 2 * np.pi
        if vy2 < 0:
            theta2 = - theta2 + 2 * np.pi

        #вычисление угла соприкосновения.
        if (y1-y2) < 0:
            phi = - np.arccos((x1-x2) / r12) + 2 * np.pi
        else:
            phi = np.arccos((x1-x2) / r12)

        # Пересчет  x-компоненты скорости первой частицы
        VX1 = v1 * np.cos(theta1 - phi) * (mass1 - K * mass2) \
        * np.cos(phi) / (mass1 + mass2)\
        + ((1 + K) * mass2 * v2 * np.cos(theta2 - phi))\
        * np.cos(phi) / (mass1 + mass2)\
        + K * v1 * np.sin(theta1 - phi) * np.cos(phi + np.pi / 2)

        # Пересчет y-компоненты скорости первой частицы
        VY1 = v1 * np.cos(theta1 - phi) * (mass1 - K * mass2) \
        * np.sin(phi) / (mass1 + mass2) \
        + ((1 + K) * mass2 * v2 * np.cos(theta2 - phi)) \
        * np.sin(phi) / (mass1 + mass2) \
        + K * v1 * np.sin(theta1 - phi) * np.sin(phi + np.pi / 2)

        # Пересчет x-компоненты скорости второй частицы
        VX2 = v2 * np.cos(theta2 - phi) * (mass2 - K * mass1) \
        * np.cos(phi) / (mass1 + mass2)\
        + ((1 + K) * mass1 * v1 * np.cos(theta1 - phi)) \
        * np.cos(phi) / (mass1 + mass2)\
        + K * v2 * np.sin(theta2 - phi) * np.cos(phi + np.pi / 2)

        # Пересчет y-компоненты скорости второй частицы
        VY2 = v2 * np.cos(theta2 - phi) * (mass2 - K * mass1) \
        * np.sin(phi) / (mass1 + mass2) \
        + ((1 + K) * mass1 * v1 * np.cos(theta1 - phi)) \
        * np.sin(phi) / (mass1 + mass2)\
        + K * v2 * np.sin(theta2 - phi) * np.sin(phi + np.pi / 2)

    else:
        # если условие столкновнеия не выполнено, то скорости частиц не пересчитываются
        VX1, VY1, VX2, VY2 = vx1, vy1, vx2, vy2

    return VX1, VY1, VX2, VY2
