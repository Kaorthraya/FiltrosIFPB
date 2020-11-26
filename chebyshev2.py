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

class f_chebyshev2:
    def __init__(self, tipo, **kwargs):
        self.a_p = kwargs['a_p']
        self.a_s = kwargs['a_s']
        self.w_s = kwargs.get('w_s', 0)
        self.w_p = kwargs.get('w_p', 0)
        self.w_s1 = kwargs.get('w_s1', 0)
        self.w_s2 = kwargs.get('w_s2', 0)
        self.w_p1 = kwargs.get('w_p1', 0)
        self.w_p2 = kwargs.get('w_p2', 0)
        self.ativo = False
        self.fcn = 0
        self.w = 0
        self.amp = 0
        self.fase = 0
        self.tipo = tipo
        self.zer = 0
        self.Rl = 0
        self.top = ''
        self.roots = 0
        self.den_norm = 0
        self.num_norm = 0
        self.Bp = 0
        self.Bs = 0
        self.eps = 0

    def ganho_bp(self, G):
        if(G != 0):
            self.G_bp = -G
            self.a_s = self.a_s + self.G_bp
            self.a_p = self.a_p + self.G_bp
            self.ativo = True
        else:
            raise ValueError("A atenuação na banda de passagem deve ser menor que o ganho na banda de passagem")

    def ordem(self, **kwargs):
        result = 0
        self.testa_par()
        if self.tipo == 'hp':
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.w_p/self.w_s))))
        if self.tipo == 'lp':
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.w_s/self.w_p))))
        if(self.tipo == 'bp'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.Bs/self.Bp))))
        if(self.tipo == 'bs'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.arccosh((np.sqrt(np.power(10, -self.a_s/10)-1))/(np.sqrt(np.power(10, -self.a_p/10)-1)))/(np.arccosh(self.Bp/self.Bs))))
        result = kwargs.get('ord', result)
        return result
        
    def fq_p(self):
        result = None
        if(self.tipo == 'hp'):
            result = self.w_s
        elif(self.tipo == 'lp'):
            result = self.w_s
        elif(self.tipo == 'bp'):
            result = np.sqrt(self.w_p1*self.w_p2)
        elif(self.tipo == 'bs'):
            result = np.sqrt(self.w_s1*self.w_s2)
        self.ws = result
        return result

    def raizes_normal(self, **kwargs):
        if ('ordem' in kwargs):
            ordem = kwargs['ordem']
        else:
            ordem = self.ordem()
        if(ordem%2 != 0):
            k = int((ordem-1)/2)
        else:
            k = int(ordem/2)
        
        self.eps = 1/np.sqrt(pow(10, (-0.1*self.a_s)) - 1)
        zeros = np.zeros(k, dtype=complex)
        alfa_k = np.zeros(ordem, dtype=complex)
        w_k = np.zeros(ordem, dtype=complex)
        D_k = np.zeros(ordem)

        for i in range(1, ordem + 1):
            D_k = (np.power(np.sinh((1/ordem)*np.arcsinh(1/self.eps)), 2)*np.power((np.sin((np.pi/(2*ordem))*(2*i-1))),2) ) \
            +(np.power(np.cosh((1/ordem)*np.arcsinh(1/self.eps)),2)*np.power(np.cos((np.pi/(2*ordem))*(2*i-1)),2))
            alfa_k[i-1] = (-1*np.sinh((1/ordem)*np.arcsinh(1/self.eps))*np.sin( (np.pi/(2*ordem))*(2*i - 1)))/D_k
            w_k[i-1] = (-1*np.cosh((1/ordem)*np.arcsinh(1/self.eps))*np.cos((np.pi/(2*ordem))*(2*i - 1)))/D_k
        self.roots = alfa_k + 1j*w_k
        
        for p in np.arange(1, k+1):
            zeros[p-1] = (1j)/(np.cos((np.pi/2)*((2*p-1)/ordem)))
        self.zer = np.append(zeros, -1*zeros)
        
        return (self.roots, self.zer)

    def transfunc(self, zeros, polos, **kwargs):
        ws = kwargs.get('ws', 0)
        w0 = kwargs.get('w0', 0)
        Bw = kwargs.get('bw', self.Bs)
        G_db = kwargs.get('G', self.G_bp)
        ordem = kwargs.get('ord', self.ordem())
        fcn = 0
        if self.tipo == 'lp':
            if(len(zeros) == 1):
                self.num_norm = np.array([1])
            else:
                self.num_norm = np.real(np.poly(zeros))
            self.den_norm = np.real(np.poly(polos))
            numm = self.num_norm
            denm = self.den_norm
            for i in range(0, len(polos) + 1):
                denm[i] = denm[i]*pow(ws, i)
            for i in range(0, len(zeros)+1):
                numm[i] = numm[i]*pow(ws, i)
            if(self.ativo == True and len(polos) == 2):
                k = 1
            else:
                k = denm[-1]/numm[-1]
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*k*numm, denm)
        if(self.tipo == 'hp'):
            self.den_norm = np.real(np.poly(polos))
            self.num_norm = np.real(np.poly(zeros))
            numm = self.num_norm[::-1]
            denm = self.den_norm[::-1]
            for i in range(0, len(polos) + 1):
                denm[i] = denm[i]*pow(ws, i)
            for i in range(0, len(zeros)+1):
                numm[i] = numm[i]*pow(ws, i)
            if(ordem%2 != 0):
                numm = np.append(numm, 0)
            k = denm[0]/numm[0]
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*k*numm, denm)
        if(self.tipo == 'bp'):
            self.den_norm = np.real(np.poly(polos))
            self.num_norm = np.real(np.poly(zeros))
            k = self.den_norm[-1]/self.num_norm[-1]
            numm, denm = signal.lp2bp(k*self.num_norm, self.den_norm, w0, Bw)
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*numm, denm)
        self.fcn = fcn
        return fcn

    def transfunc2(self, t):
        ordem = self.ordem()
        polos2o = np.zeros(2)
        zeros2o = np.zeros(2)
        zerosNo = np.roots(self.fcn.num)
        polosNo = np.roots(self.fcn.den)
        modZ = np.abs(zerosNo)
        modP = np.abs(polosNo)
        zerosNo = zerosNo[np.argsort(modZ)]
        polosNo = polosNo[np.argsort(modP)]
        tf2o = []
        if(t.lower() == 'notch' and (self.tipo == 'lp' or self.tipo == 'hp')):
            if(ordem%2 == 0):
                for k in range(0, int(ordem/2)):
                    zeros2o = [zerosNo[2*k], zerosNo[2*k + 1]]
                    polos2o = [polosNo[2*k], polosNo[2*k + 1]]
                    print(zeros2o)
                    print(polos2o)
                    num = np.poly((1j)*np.imag(zeros2o))
                    den = np.poly(polos2o)
                    tf = signal.TransferFunction(num, den)
                    tf2o.append(tf)
            else:
                for k in range(0, int(ordem/2)+1):
                    if(k != int(ordem/2)):
                        zeros2o = [zerosNo[2*k], zerosNo[2*k + 1]]
                        polos2o = [polosNo[2*k], polosNo[2*k + 1]]
                        num = np.poly((1j)*np.imag(zeros2o))
                        den = np.poly(polos2o)
                        tf = signal.TransferFunction(num, den)
                        tf2o.append(tf)
                    else:
                        zero1o = 1
                        polo1o = polosNo[-1]
                        num = np.poly((1j)*np.imag(zero1o))
                        den = np.poly(polo1o)
                        tf = signal.TransferFunction(num, den)
                        tf2o.append(tf)
        self.tf2od = tf2o
        return tf2o

    def graphpoints(self, min_f, max_f, points):
        self.w, self.amp, self.fase = self.fcn.bode(w=np.arange(min_f, max_f, (max_f - min_f)/points))

    def plot_bode(self, **kwargs):
        
        #FIGURA 1 - AMPLITUDE (dB)
        fig = plt.figure(1)
        plt.semilogx(self.w, self.amp)
        ax = plt.gca()
        ax.set(xlabel="Frequência (rad/s)", ylabel="Amplitude em dB",
               title="Resposta em amplitude (CHEBYSHEV TIPO 2: n = %d)" %self.ordem()) # configuração de plot label
        ax.margins(x=0)
        ax.margins(y=0.05)  # margem y
        ax.grid(True, which="both")  # grid

        if(self.tipo == 'hp'  or self.tipo == 'lp'):
            bp = ax.scatter(self.w_p, self.a_p - self.G_bp)  # ponto de projeto de passagem
            br = ax.scatter(self.w_s, self.a_s - self.G_bp)  # ponto de projeto de rejeição
            ax.legend((bp, br), ("P. Projeto (passagem)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        if(self.tipo == 'bp' or self.tipo == 'bs'):
            bp1 = ax.scatter(self.w_p1, self.a_p)
            bp2 = ax.scatter(self.w_p2, self.a_p)  # ponto de projeto de passagem
            bs1 = ax.scatter(self.w_s1, self.a_s)  # ponto de projeto de rejeição
            bs2 = ax.scatter(self.w_s2, self.a_s)
            ax.legend((bp1, bp2, bs1, bs2), ("P. Projeto (passagem)", "P. Projeto (passagem)", "P. Projeto (rejeição)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        
        #FIGURA 2 - FASE
        fig2 = plt.figure(2)
        plt.semilogx(self.w, self.fase)
        ax2 = plt.gca()
        ax2.set(xlabel="Frequência (rad/s)", ylabel="Fase em graus",
               title="Resposta em fase (CHEBYSHEV TIPO 2: n = %d)" %self.ordem()) # configuração de plot label
        ax2.margins(x=0)
        ax2.margins(y=0.05)  # margem y
        ax2.grid(True, which="both")  # grid

        return fig, fig2

    def testa_par(self):
        if (self.tipo == 'hp') and ((self.w_s >= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um HP")
        if (self.tipo == 'lp') and ((self.w_s <= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um LP")
