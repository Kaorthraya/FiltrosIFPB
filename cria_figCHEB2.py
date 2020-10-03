import numpy as np
from matplotlib import pyplot as plt
from chebyshev2 import *

at = np.arange(1,11,1, dtype=int)
print(at)
for i in range(0, 10):
    filtro = f_chebyshev2('lp', a_p = -1, a_s = -40, w_p = 0.2, w_s = 1)
    print(at[i])
    n = filtro.ordem(ord = at[i])
    fc1 = filtro.fq_corte()
    p, z = filtro.raizes_normal(ordem = n)
    fc = filtro.transfunc(z, p, wc = fc1, ord = at[i])
    filtro.graphpoints(fc, 10, 0.001, 10000)
    fig = plt.figure(1)
    plt.semilogx(filtro.w, filtro.amp, label = 'Ordem: %d' %n)
    ax = plt.gca()
    ax.set(xlabel="Frequência (rad/s)", ylabel="Amplitude em dB",
               title="Resposta em amplitude (CHEBYSHEV TIPO 2)")  # configuração de plot label
    ax.margins(x=0)
    ax.margins(y=0.05)  # margem y
    ax.grid(True, which="both")  # grid
    ax.set_ylim(-80, 1)
plt.savefig('CHEB2.png', dpi = 1000)
plt.legend()
plt.show()

