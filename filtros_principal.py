# -*- coding: utf-8 -*-

from cl_fc import f_btrwrth, eng_string

import numpy as np
from matplotlib import pyplot as plt

print("Tipos de respostas calculadas: Highpass [hp] ou Lowpass [lp]")
resposta = input("Qual o tipo do filtro que deseja: ")
print("Pontos de projeto: ")
w_s = float(input("Digite a frequência de rejeição: "))
a_s = float(input("Atenuação na frequência de rejeição: "))
w_p = float(input("Digite a frequência de passagem: "))
a_p = float(input("Atenuação na frequência de passagem: "))

filtro = f_btrwrth(resposta, a_p = a_p, w_p = w_p, a_s = a_s, w_s = w_s)
n = filtro.ordem()
fc1 = filtro.fq_corte()
raiz1 = filtro.raizes_normal()
tf1 = filtro.transfunc(raiz1,fc1)

choice = input("Deseja fazer o escalonamento em frequência? [y/n]: ")
if(choice == 'y'):
    R = float(input("Valor da carga (ohms): "))
    print("Qual topologia deseja usar: ")
    print("[CLC] -> Começa com capacitor\n[LCL] -> Começa com indutor")
    top = input("Topologia: ")
    valor, nome = filtro.elements(fc1, top, R)
    print("Elementos escalonados:\n")
    for i in range(0, n):
        print('%s= %s' %(nome[i], eng_string(valor[i], si=True)))
    print("\n")

print("Frequência de corte: %f" %fc1)

filtro.exporttxt('saida.txt')

graf1 = filtro.plot_bode(tf1, min_f=1,max_f=10000, points=1000)
graf1.show()

graf1.savefig('fig1.png', dpi = 600)
plt.show()

