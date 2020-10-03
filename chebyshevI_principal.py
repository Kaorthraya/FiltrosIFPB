# -*- coding: utf-8 -*-

from chebyshev2 import *

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

print('Cálculo de filtros Butterworth\n')
print('Opções de resposta:\n[lp] - lowpass\n[hp] - highpass\n[bp] - bandpass\n[bs] - bandstop\n')
resposta = input('Digite a resposta do filtro: ')
if(resposta.lower() == 'lp' or resposta.lower() == 'hp'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp = float(input('Valor de frequência (passagem): '))
    Ws = float(input('Valor de frequência (rejeição): '))
    filtro = f_chebyshev(resposta.lower(),a_p = ap_dB,
                                        a_s = as_dB,
                                        w_p = Wp,
                                        w_s = Ws)
elif(resposta.lower() == 'bp' or resposta.lower() == 'bs'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp1 = float(input('Valor inferior da banda de passagem: '))
    Wp2 = float(input('Valor superior da banda de passagem: '))
    Ws1 = float(input('Valor inferior da banda de rejeição: '))
    Ws2 = float(input('Valor superior da banda de rejeição: '))
    filtro = f_chebyshev(resposta.lower(),a_p = ap_dB,
                                        a_s = as_dB,
                                        w_s1 = Ws1,
                                        w_s2 = Ws2,
                                        w_p1 = Wp1,
                                        w_p2 = Wp2)

n = filtro.ordem()
print('Ordem: %d' %n)
fc1 = filtro.fq_corte()
polos, zeros = filtro.raizes_normal()
tf = filtro.transfunc(zeros, polos, wc = fc1)
print(fc1)
print(tf)
filtro.graphpoints(tf, 100e3, 1, 100e3)
filtro.plot_bode()
plt.show()

#if(resposta == 'lp' or resposta == 'hp'):
#    tf1 = filtro.transfunc(raiz1, wc = fc1)
#else:
#    tf1 = filtro.transfunc(raiz1, w0 = fc1)

#graf1 = filtro.plot_bode(tf1, min_f=1,max_f=1e6, points=100e3)
#graf1.canvas.set_window_title('TF')

#plt.show()




