import numpy as np
from matplotlib import pyplot as plt
from butterworth import *
from chebyshev import *
from chebyshev2 import *

filtro_bt = f_btrwrth('bp', w_p1 = 2000 , w_p2 = 5000, w_s1 = 1000, w_s2 = 10000, a_p = -1, a_s = -40)
filtro_chb1 = f_chebyshev1('bp', w_p1 = 2000 , w_p2 = 5000, w_s1 = 1000, w_s2 = 10000, a_p = -1, a_s = -40)
filtro_chb2 = f_chebyshev2('bp', w_p1 = 2000 , w_p2 = 5000, w_s1 = 1000, w_s2 = 10000, a_p = -1, a_s = -40)

n_bt = filtro_bt.ordem()
n_chb1 = filtro_chb1.ordem()
n_chb2 = filtro_chb2.ordem()

fc_bt = filtro_bt.fq_corte()
fc_chb1 = filtro_chb1.fq_p()
fc_chb2 = filtro_chb2.fq_p()

r_bt = filtro_bt.raizes_normal()
p_chb1 = filtro_chb1.raizes_normal()
p_chb2, z_chb2 = filtro_chb2.raizes_normal()

tf_bt = filtro_bt.transfunc(r_bt, w0 = fc_bt)
tf_chb1 = filtro_chb1.transfunc(p_chb1, w0 = fc_chb1)
tf_chb2 = filtro_chb2.transfunc(z_chb2, p_chb2, w0 = fc_chb2)

filtro_bt.graphpoints(tf_bt, 500e3, 10, 500e3)
filtro_chb1.graphpoints(500e3, 10, 500e3)
filtro_chb2.graphpoints(500e3, 10, 500e3)

fig1 = plt.figure(1)
plt.semilogx(filtro_bt.w, filtro_bt.amp, label = 'Butterworth: n = %d' %n_bt)
plt.semilogx(filtro_chb1.w, filtro_chb1.amp, label = 'Chebyshev I: n = %d' %n_chb1)
plt.semilogx(filtro_chb2.w, filtro_chb2.amp, label = 'Chebyshev II: n = %d' %n_chb2)
ax = plt.gca()
ax.set(xlabel="Frequência (rad/s)", ylabel="Amplitude em dB",
    title="Resposta em amplitude" ) # configuração de plot label
ax.margins(x=0)
ax.margins(y=0.05)  # margem y
ax.grid(True, which="both")  # grid
ax.set_ylim(-80, 5)
ax.legend()

fig1 = plt.figure(2)
plt.semilogx(filtro_bt.w, filtro_bt.fase, label = 'Butterworth: n = %d' %n_bt)
plt.semilogx(filtro_chb1.w, filtro_chb1.fasedeg, label = 'Chebyshev I: n = %d' %n_chb1)
plt.semilogx(filtro_chb2.w, filtro_chb2.fase, label = 'Chebyshev II: n = %d' %n_chb2)
ax = plt.gca()
ax.set(xlabel="Frequência (rad/s)", ylabel="Fase em graus",
    title="Resposta em fase" ) # configuração de plot label
ax.margins(x=0)
ax.margins(y=0.05)  # margem y
ax.grid(True, which="both")  # grid
ax.legend()


plt.show()
