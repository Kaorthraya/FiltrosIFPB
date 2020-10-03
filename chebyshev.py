# -*- coding: utf-8 -*-
"""cl_fc.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1vOo2SJUOSCdwN6rL7sUBGzzABnTlzVeM
"""

import numpy as np
import math
from scipy import signal
from matplotlib import pyplot as plt
from matplotlib.ticker import (AutoMinorLocator)

class f_chebyshev1:
    def __init__(self, tipo, **kwargs):
        self.a_p = kwargs['a_p']
        self.a_s = kwargs['a_s']
        self.w_s = kwargs.get('w_s', 0)
        self.w_p = kwargs.get('w_p', 0)
        self.w_s1 = kwargs.get('w_s1', 0)
        self.w_s2 = kwargs.get('w_s2', 0)
        self.w_p1 = kwargs.get('w_p1', 0)
        self.w_p2 = kwargs.get('w_p2', 0)
        self.tipo = tipo
        self.gpdelay = 0
        self.fasedeg = 0
        self.w = 0
        self.amp = 0
        self.fase = 0
        self.wc = 0.0
        self.Rl = 0
        self.top = ''
        self.roots = 0
        self.den_norm = 0
        self.Bp = 0
        self.Bs = 0
        self.eps = 0

    def ordem(self, **kwargs):
        result = 0
        self.testa_par()
        if self.tipo == 'hp':
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.w_p/self.w_s))))
        if self.tipo == 'lp':
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.w_s/self.w_p))))
        if(self.tipo == 'BP'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.Bs/self.Bp))))
        if(self.tipo == 'bs'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.Bp/self.Bs))))
        result = kwargs.get('ord', result)
        return result
        
    def fq_corte(self):
        result = None
        ordem = self.ordem()
        if(self.tipo == 'hp'):
            result = self.w_p
        elif(self.tipo == 'lp'):
            result = self.w_p
        elif(self.tipo == 'BP'):
            result = np.sqrt(self.w_p1*self.w_p2)
        elif(self.tipo == 'bs'):
            result = np.sqrt(self.w_s1*self.w_s2)
        self.wc = result
        return result

    def raizes_normal(self, **kwargs):
        if ('ordem' in kwargs):
            ordem = kwargs['ordem']
        else:
            ordem = self.ordem()
        self.eps = np.sqrt(np.power(10, -self.a_p/10) - 1)
        alfa_k = np.zeros(ordem, dtype=complex)
        w_k = np.zeros(ordem, dtype=complex)
        for i in range(1, ordem + 1):
            alfa_k[i-1] = -1*np.sinh((1/ordem)*np.arcsinh(1/self.eps))*np.sin( (np.pi/(2*ordem))*(2*i - 1))
            w_k[i-1] = 1*np.cosh((1/ordem)*np.arcsinh(1/self.eps))*np.cos((np.pi/(2*ordem))*(2*i - 1))
        self.roots = alfa_k + 1j*w_k
        return self.roots

    def transfunc(self, polos, **kwargs):
        wc = kwargs.get('wc', 0)
        w0 = kwargs.get('w0', 0)
        Bw = kwargs.get('bw', self.Bp)
        ordem = kwargs.get('ord', self.ordem())
        fcn = 0
        if self.tipo == 'lp':
            self.den_norm = np.real(np.poly(polos))
            denm = np.zeros(len(self.den_norm))
            for i in range(0, len(polos) + 1):
                denm[i] = self.den_norm[i] * np.power(wc, i)
            if(ordem%2 == 0):
                num = denm[-1]*(1/np.sqrt(1+np.power(self.eps,2)))
            else: 
                num = denm[-1]
            fcn = signal.TransferFunction(num, denm)
        if self.tipo == 'hp':
            self.den_norm = np.real(np.poly(polos))
            denm = np.zeros(len(polos) + 1)
            for i in range(0, len(polos) + 1):
                denm[i] = self.den_norm[len(polos) - i] * np.power(wc, i)
            num = np.zeros(len(polos) + 1)
            if(self.ordem()%2 == 0):
                num = denm[-1]*(1/np.sqrt(1+np.power(self.eps,2)))
            else: 
                num[0] = denm[0]
            fcn = signal.TransferFunction(num, denm)
        if(self.tipo == 'BP'):
            self.den_norm = np.real(np.poly(polos))
            [num, den] = signal.lp2bp(self.den_norm[-1], self.den_norm, w0, Bw)
            if(self.ordem()%2 == 0):
                num[0] = num[0]*(1/np.sqrt(1+np.power(self.eps,2)))
            else: 
                num[0] = num[0]
            fcn = signal.TransferFunction(num, den)
        if(self.tipo == 'bs'):
            self.den_norm = np.real(np.poly(polos))
            [num, den] = signal.lp2bs(self.den_norm[-1], self.den_norm, w0, Bw)
            if(self.ordem()%2 == 0):
                num = num[::]*(1/np.sqrt(1+np.power(self.eps,2)))
            else: 
                num[-1] = num[-1]
            fcn = signal.TransferFunction(num, den)

        return fcn

    def graphpoints(self, fcn, max_f, min_f, points):
        self.w, self.amp, self.fasedeg = fcn.bode(w=np.arange(min_f, max_f, (max_f - min_f)/points))
        fase_aj = np.zeros(len(self.fasedeg))
        for i in range(0, len(self.fasedeg)):
            fase_aj[i] = self.fasedeg[i]*(np.pi/180.0) #AJUSTA A FASE PARA RADIANOS
        self.fase = fase_aj

    def gp_delay(self):
        aux = ((-1)*np.diff(self.fase)/np.diff(self.w))
        self.gpdelay = np.append(aux, 0)        #ADICIONA UM ZERO NO FIM PARA CASAR OS TAMANHOS

    def plot_bode(self, fcn, **kwargs):

        fig, ax = plt.subplots()  # cria os plots
        ax.semilogx(self.w, self.amp)  # gráfico do tipo semilog
        ax.set(xlabel="Frequência (rad/s)", ylabel="Amplitude em dB",
               title="Resposta em amplitude (BTRWRTH n = %d)" %self.ordem()) # configuração de plot label
        ax.margins(x=0)
        ax.margins(y=0.05)  # margem y
        ax.grid(True, which="both")  # grid
        if(self.tipo == 'hp'  or self.tipo == 'lp'):
            bp = ax.scatter(self.w_p, self.a_p)  # ponto de projeto de passagem
            br = ax.scatter(self.w_s, self.a_s)  # ponto de projeto de rejeição
            ax.legend((bp, br), ("P. Projeto (passagem)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        if(self.tipo == 'BP' or self.tipo == 'bs'):
            bp1 = ax.scatter(self.w_p1, self.a_p)
            bp2 = ax.scatter(self.w_p2, self.a_p)  # ponto de projeto de passagem
            bs1 = ax.scatter(self.w_s1, self.a_s)  # ponto de projeto de rejeição
            bs2 = ax.scatter(self.w_s2, self.a_s)
            ax.legend((bp1, bp2, bs1, bs2), ("P. Projeto (passagem)", "P. Projeto (passagem)", "P. Projeto (rejeição)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        return fig

    def testa_par(self):
        if (self.tipo == 'hp') and ((self.w_s >= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um HP")
        if (self.tipo == 'lp') and ((self.w_s <= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um LP")
