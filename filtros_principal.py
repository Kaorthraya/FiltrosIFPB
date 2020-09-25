# -*- coding: utf-8 -*-

#from chebyshev import *
from butterworth import *

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
    filtro = f_btrwrth(resposta.lower(),a_p = ap_dB,
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
    filtro = f_btrwrth(resposta.lower(),a_p = ap_dB,
                                        a_s = as_dB,
                                        w_s1 = Ws1,
                                        w_s2 = Ws2,
                                        w_p1 = Wp1,
                                        w_p2 = Wp2)

n = filtro.ordem()
fc1 = filtro.fq_corte()
raiz1 = filtro.raizes_normal()

print('Ordem: %d' %n)

if(resposta == 'lp' or resposta == 'hp'):
    tf1 = filtro.transfunc(raiz1, wc = fc1)
else:
    tf1 = filtro.transfunc(raiz1, w0 = fc1)

graf1 = filtro.plot_bode(tf1, min_f=1,max_f=1e6, points=100e3)
graf1.canvas.set_window_title('TF')


while(filtro.criterio == -1):
    print('Algum critério não foi atendido\n')
    print('Modificação: aumento ou diminuição de banda')
    ch1 = input('Deseja realizar otimização na largura de banda? [s/n]:')
    flag_opt = 0
    if(ch1.lower == 's'):
        pass
    else:
        print('Otimização não realizada. Parâmetros originais preservados\n')
        break
    pct = float(input('Digite o valor em porcentagem do aumento/diminuição de banda: '))
    tf, n_bw = filtro.bw_increase(pct)
    graf1 = filtro.plot_bode(tf, min_f=1,max_f=1e6, points=100e3)
    graf1.canvas.set_window_title('Banda aumentada em %2.1f %%' %pct)
    print('Novo valor da banda: %2.1f' %n_bw)
    flag_opt = 1
    plt.show()

if(resposta.lower() == 'bp' or resposta.lower() == 'bs'):
    if(flag_opt == 1):
        print('BW a ser usado para calcular os elementos:\n[1] - Original\n[2] - Otimizado\n')
        ch = int(input('Digite sua escolha: '))
        if(ch == 1):
            n_bw = Wp2 - Wp1
        elif(ch == 2):
            n_bw = n_bw
    else:
        n_bw = Wp2 - Wp1

    print('Escolha a topologia: [S] -> Começa com LC série\n[P] -> Começa com LC paralelo\n')
    top = input('Digite sua escolha: ')
    R = float(input('Digite o valor da impedância (Rs = RL): '))
    a, b = filtro.elements(fc1, top.lower(), R, bw = n_bw)
    for i in range(0, 2*n):
        print('     %s = %s' %(b[i], eng_string(a[i], si=True)))    

graf1.savefig('fig1.png', dpi = 600)



