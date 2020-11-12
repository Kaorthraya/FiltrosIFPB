# -*- coding: utf-8 -*-

from chebyshev2 import *
from printft import printTF

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

print('Cálculo de filtros Chebyshev 2\n')
print('Opções de resposta:\n[lp] - lowpass\n[hp] - highpass\n[bp] - bandpass\n[bs] - bandstop\n')
resposta = input('Digite a resposta do filtro: ')
if(resposta.lower() == 'lp' or resposta.lower() == 'hp'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp = float(input('Valor de frequência (passagem): '))
    Ws = float(input('Valor de frequência (rejeição): '))
    filtro = f_chebyshev2(resposta.lower(),a_p = ap_dB,
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
    filtro = f_chebyshev2(resposta.lower(),a_p = ap_dB,
                                        a_s = as_dB,
                                        w_s1 = Ws1,
                                        w_s2 = Ws2,
                                        w_p1 = Wp1,
                                        w_p2 = Wp2)

if(ap_dB >= 0):
    Ganho = float(input('Ganho na banda de passagem (dB): '))
    filtro.ganho_bp(Ganho)
else:
    filtro.ganho_bp(0)

n = filtro.ordem()
print('Ordem: %d' %n)
fc1 = filtro.fq_p()
polos, zeros = filtro.raizes_normal()

if(resposta == 'lp' or resposta == 'hp'):
    tf = filtro.transfunc(zeros, polos, ws = fc1)
    printTF(tf.num, tf.den)
    filtro.graphpoints(1,100e3,100e3)
    filtro.plot_bode()
else:
    tf1 = filtro.transfunc(zeros, polos, w0 = fc1)
    filtro.graphpoints(1, 100e3, 100e3)
    filtro.plot_bode()

if(n >= 2):
    tfs2o = filtro.transfunc2('notch')
    for i in range(0, len(tfs2o)):
        printTF(tfs2o[i].num, tfs2o[i].den)

plt.show()




