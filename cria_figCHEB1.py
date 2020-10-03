import numpy as np
from matplotlib import pyplot as plt
from chebyshev import *

#PARA ALTERAR PARA FASE, ALTERAR SEMILOGX PARA A VARIAVEL FILTRO.GPDELAY

at = np.arange(1,11,1, dtype=int)

for i in range(0, 10):
    filtro = f_chebyshev('lp', a_p = -1, a_s = -40, w_p = 1, w_s = 2)
    n = filtro.ordem(ord = at[i])
    fc1 = filtro.fq_corte()
    p = filtro.raizes_normal(ordem = n)
    fc = filtro.transfunc(p, wc = fc1, ord = at[i])
    filtro.graphpoints(fc, 10, 0.1, 100e3)
    filtro.gp_delay()
    fig = plt.figure(1)
    plt.semilogx(filtro.w, filtro.gpdelay, label = 'Ordem: %d' %n)
    ax = plt.gca()
    ax.set(xlabel="Frequência (rad/s)", ylabel="Delay de grupo (s)",
               title="Delay de grupo (CHEBYSHEV TIPO 1)")  # configuração de plot label
    ax.margins(x=0)
    ax.margins(y=0.05)  # margem y
    ax.grid(True, which="both")  # grid
plt.savefig('CHEB2.png', dpi = 1000)
plt.legend()
plt.show()

