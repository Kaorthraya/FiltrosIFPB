# -*- coding: utf-8 -*-

from chebyshev import *
from printft import printTF

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

print('Cálculo de filtros Chebyshev I\n')
print('Opções de resposta:\n[lp] - lowpass\n[hp] - highpass\n[bp] - bandpass\n[bs] - bandstop\n')
resposta = input('Digite a resposta do filtro: ')
if(resposta.lower() == 'lp' or resposta.lower() == 'hp'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp = float(input('Valor de frequência (passagem): '))
    Ws = float(input('Valor de frequência (rejeição): '))
    filtro = f_chebyshev1(resposta.lower(),a_p = ap_dB,
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
    filtro = f_chebyshev1(resposta.lower(),a_p = ap_dB,
                                        a_s = as_dB,
                                        w_s1 = Ws1,
                                        w_s2 = Ws2,
                                        w_p1 = Wp1,
                                        w_p2 = Wp2)

tipo = input('O filtro é ativo [A] ou passivo [P]: ')
if(tipo.lower() == 'a'):
    G = float(input('Ganho na banda de passagem (dB): '))
    filtro.ganho_bp(G)
else:
    filtro.ganho_bp(0)

n = filtro.ordem()
print('Ordem: %d' %n)
fc1 = filtro.fq_p()
polos = filtro.raizes_normal()


if(resposta.lower() == 'lp' or resposta.lower() == 'hp'):
    tf = filtro.transfunc(polos, wp = fc1)
    filtro.graphpoints(1,1e6,1e6)
    print('Frequência de corte: %f' %filtro.wp)
    filtro.plot_bode()
else:
    tf = filtro.transfunc(polos, w0 = fc1)
    filtro.graphpoints(100e3, 1, 100e3)
    print('Frequência de corte: %f' %filtro.wp)
    filtro.plot_bode()

print("\nFunção de transferência do filtro: ")
printTF(tf.num, tf.den)

tfs2o = filtro.transfunc2('tow')

if(n%2 == 0):
    for i in range(0, int(n/2)):
        printTF(tfs2o[i].num, tfs2o[i].den)
        Q = np.sqrt(tfs2o[i].den[-1])/(tfs2o[i].den[-2])
        print('Q = %2.2f' %Q)
else:
    for i in range(0, int(n/2)+1):
        printTF(tfs2o[i].num, tfs2o[i].den)
        if(i != int(n/2)):
            Q = np.sqrt(tfs2o[i].den[-1])/(tfs2o[i].den[-2])
            print('Q = %2.2f' %Q)

plt.show()




