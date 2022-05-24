# Autor: Douglas Mendes
# Data: 24/maio/2022
# Variação das propriedades termodinâmicas no ar em uma onda de choque normal
# em função da velocidade Mach (gamma = 1.2)


import numpy as np
import math
import matplotlib.pyplot as plt

# Definição das Constantes
gamma = 1.2
R = 8.314462618     # J * kg−1 * K−1

def calc_m2(m1):
    m1_squared = m1**2
    m2 = ((1 + (((gamma-1)/2)*m1_squared)) / (gamma*m1_squared-((gamma-1)/2)))
    m2 = math.sqrt(m2)
    return m2

def calc_razao_rho(m1):
    m1_squared = m1**2
    razao_rho = ((gamma+1)*m1_squared) / (2+((gamma-1)*m1_squared))
    return razao_rho 

def calc_razao_p(m1):
    m1_squared = m1**2
    razao_p = (1+((2*gamma)/(gamma+1))*(m1_squared-1))
    return razao_p

def calc_razao_T(m1):
    m1_squared = m1**2
    razao_p = calc_razao_p(m1) 
    inv_razao_rho = (2+((gamma-1)*m1_squared)) / ((gamma+1)*m1_squared)
    razao_T = razao_p * inv_razao_rho
    return razao_T

def calc_razao_p01_p02(m1):
    dif_entropia = calc_razao_p(m1)  
    razao_p01_p02 = math.exp(((-1)*dif_entropia)/R)
    return razao_p01_p02

mach_1 = np.arange(1.0, 10.0, 0.01)

# Cálculo de Mach 2:
mach_2 = np.array([calc_m2(xi) for xi in mach_1])

# Cálculo da razão p0,1/p0,2:
razao_p01_p02 = np.array([calc_razao_p01_p02(xi) for xi in mach_1])

# Cálculo da razão p1/p2:
razao_p1_p2 = np.array([calc_razao_p(xi) for xi in mach_1])

# Cálculo da razão rho_1/rho_2:
razao_rho1_rho2 = np.array([calc_razao_rho(xi) for xi in mach_1])

# Cálculo da razão de Temperaturas:
razao_T1_T2 = np.array([calc_razao_T(xi) for xi in mach_1])

# -----------------------------
# Lado esquerdo:
fig, plot_esq = plt.subplots()


plot_esq.set_xlabel('M1')
plot_esq.set_ylabel('M2 e p0,2/p0,1', color='tab:red')
plot_esq.set_ylim(0.0,1.0)

color_m2 = 'tab:red'
plot_esq.plot(mach_1, mach_2, color=color_m2, label="M2")
plot_esq.tick_params(axis='y', labelcolor=color_m2)

color_p = 'tab:green'
plot_esq.plot(mach_1, razao_p01_p02, color=color_p, label="p0,1/p0,2")
plot_esq.tick_params(axis='y', labelcolor=color_p)
plot_esq.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=3, fancybox=True)

# -----------------------------
# Lado direito:

plot_dir = plot_esq.twinx()  # inicia um segundo eixo compartilhando o mesmo eixo-x
plot_dir.set_ylabel('T2/T1, p2/p1 e rho2/rho1', color='tab:blue')  # o rótulo de x já foi criado em plot_esq
plot_dir.set_ylim(0.0,20.0)

color_raz_p = 'tab:blue'
plot_dir.plot(mach_1, razao_p1_p2, color=color_raz_p, label="p1/p2")
plot_dir.tick_params(axis='y', labelcolor=color_raz_p)

color_raz_T = 'tab:orange'
plot_dir.plot(mach_1, razao_T1_T2, color=color_raz_T, label="T1/T2")
plot_dir.tick_params(axis='y', labelcolor=color_raz_T)

color_raz_rho = 'tab:purple'
plot_dir.plot(mach_1, razao_rho1_rho2, color=color_raz_rho, label="rho1/rho2")
plot_dir.tick_params(axis='y', labelcolor=color_raz_rho)

plot_dir.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),
          fancybox=True, ncol=5)

# -----------------------------
# Exibição final:
fig.tight_layout()  # caso contrário tira-se um pouco do rótulo y direito
fig.savefig('gamma_1.2.png')
plt.show()
