# -*- coding: utf-8 -*-

from chebyshev import *
#from butterworth import *

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

print('Cálculo de filtros Butterworth\n')
print('Opções de resposta:\n[lp] - lowpass\n[hp] - highpass\n[bp] - bandpass\n[bs] - bandstop\n')
resposta = input('Digite a resposta do filtro: ')
if(resposta == 'lp' or resposta == 'hp'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp = float(input('Valor de frequência (passagem): '))
    Ws = float(input('Valor de frequência (rejeição): '))
    filtro = f_chebyshev(resposta,a_p = ap_dB,
                                a_s = as_dB,
                                w_p = Wp,
                                w_s = Ws)
elif(resposta == 'bp' or resposta == 'bs'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp1 = float(input('Valor inferior da banda de passagem: '))
    Wp2 = float(input('Valor superior da banda de passagem: '))
    Ws1 = float(input('Valor inferior da banda de rejeição: '))
    Ws2 = float(input('Valor superior da banda de rejeição: '))
    filtro = f_chebyshev(resposta,a_p = ap_dB,
                                a_s = as_dB,
                                w_s1 = Ws1,
                                w_s2 = Ws2,
                                w_p1 = Wp1,
                                w_p2 = Wp2)

n = filtro.ordem()
fc1 = filtro.fq_corte()
raiz1 = filtro.raizes_normal()

if(resposta == 'lp' or resposta == 'hp'):
    tf1 = filtro.transfunc(raiz1, wc = fc1)
else:
    tf1 = filtro.transfunc(raiz1, w0 = fc1)

print(n)

graf1 = filtro.plot_bode(tf1, min_f=1,max_f=1e6, points=100e3)
plt.show()
graf1.savefig('fig1.png', dpi = 600)



