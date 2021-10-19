# Лабораторная работа №6. Задание 3
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from collisions import collision

radius = 0.5  # Радиус шариков
radius_c = 1.5
q = -10**(-7)
mass = 10  # Масса шариков

# Границы
Lx = 20
Ly = 20

K = 1 # Коэффициент столкновений между шариками

T = 10 # Общее время анимации
n = 500 # Количество итераций / кадров
tau = np.linspace(0,T,n) # Массив для одного шага анимации
dT = T / n # Время одного шага итерации

N = 8 # Количество чатсиц
p = np.zeros((N,4)) # Массивы для текущих значений положений и скоростей частиц

# Массивы для записи итоговых координат на каждой итерации для итоговой анимации
x = np.zeros((N,n))
y = np.zeros((N,n))

# Массивы для записи х, y, vx, vy для каждой частицы
p[0,0], p[0,1], p[0,2], p[0,3] = -15, 5, 6.5, 1.5
p[1,0], p[1,1], p[1,2], p[1,3] = -15, 1, 5, 1
p[2,0], p[2,1], p[2,2], p[2,3] = -15, -10, 1, 5

x[0,0], y[0,0] = p[0,0], p[0,1]
x[1,0], y[1,0] = p[1,0], p[1,1]
x[2,0], y[2,0] = p[2,0], p[2,1]


qc = 10
xc = 1
yc = 1
mass_c = 1000
const =  8.98755 * 10**9
g = 0.00  # Ускорение свободного падения

def circle_func(x_centre_point,
                y_centre_point,
                R):
    x = np.zeros(30)
    y = np.zeros(30)
    for i in range(0, 30, 1):
        alpha = np.linspace(0, 2*np.pi, 30)
        x[i] = x_centre_point + R*np.cos(alpha[i])
        y[i] = y_centre_point + R*np.sin(alpha[i])

    return x, y

def move_func(s, t):
    x, v_x, y, v_y = s

    dxdt = v_x
    dv_xdt = const*q*qc/mass*(x-xc)/((x-xc)**2+(y-yc)**2)**1.5

    dydt = v_y
    dv_ydt = const*q*qc/mass*(yc-yc)/((x-xc)**2+(y-yc)**2)**1.5

    return dxdt, dv_xdt, dydt, dv_ydt


for k in range(n-1):  # Цикл перебора шагов временеи анимации
    t = [tau[k],tau[k+1]]

    for m in range(N):  # Цикл перебора частиц
        s0 = p[m,0], p[m,2], p[m,1], p[m,3]
        sol = odeint(move_func, s0, t)

        # Перезаписываем положения частиц
        p[m,0] = sol[1,0]
        p[m,2] = sol[1,1]
        p[m,1] = sol[1,2]
        p[m,3] = sol[1,3]

        # Заноим новые положения в итоговый массив для анимации
        x[m,k+1], y[m,k+1] = p[m,0], p[m,1]

        # Проверка условий столкновения с первой перегородкой
        res = collision(p[m,0],p[m,1],p[m,2],p[m,3],xc,yc,0,0,radius,radius_c,mass,mass_c,K)
        p[m,2], p[m,3] = res[0], res[1] # Пересчет скоростей


    # Циклы перебора частиц для столкновений друг с другом
    for i in range(N): # Базовая частица
        x1, y1, vx1, vy1 = p[i,0], p[i,1], p[i,2], p[i,3] # Запись текущих координат базовой частицы
        x10, y10 = x[i,k], y[i,k] # Запись координат предыдущего шага базовой частицы

        for j in range(i+1,N): # Запись текущих координат остальных частиц
            x2, y2, vx2, vy2 = p[j,0], p[j,1], p[j,2], p[j,3] # Запись текущих
            x20, y20 = x[j,k], y[j,k] # Запись координат предыдущего шага

            # Проверка условий столкновения
            r1 = np.sqrt((x1-x2)**2+(y1-y2)**2)
            r0 = np.sqrt((x10-x20)**2+(y10-y20)**2)
            if  r1 <= radius*2 and r0 > 2*radius:
                res = collision(x1,y1,vx1,vy1,x2,y2,vx2,vy2,radius,radius,mass,mass,K)

                # Перезаписывание условий, в случае столкновения
                p[i,2], p[i,3] = res[0], res[1]
                p[j,2], p[j,3] = res[2], res[3]

# Графический вывод
fig = plt.figure()

ball1, = plt.plot([], [], 'o', color='r', ms=1)
ball2, = plt.plot([], [], 'o', color='r', ms=1)
ball3, = plt.plot([], [], 'o', color='r', ms=1)

plt.plot([xc], [yc], 'o', ms=20)

def animate(i):
    ball1.set_data(circle_func(x[0, i], y[0, i], radius))
    ball2.set_data(circle_func(x[1, i], y[1, i], radius))
    ball3.set_data(circle_func(x[2, i], y[2, i], radius))

ani = FuncAnimation(fig, animate, frames=n, interval=1)

plt.axis('equal')
plt.xlim(-Lx, Lx)
plt.ylim(-Ly, Ly)
plt.show()
