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
import warnings

class f_btrwrth:
    def __init__(self, tipo, **kwargs):
        self.a_p = kwargs['a_p']
        self.a_s = kwargs['a_s']
        self.w_s = kwargs.get('w_s', 0)
        self.w_p = kwargs.get('w_p', 0)
        self.w_s1 = kwargs.get('w_s1', 0)
        self.w_s2 = kwargs.get('w_s2', 0)
        self.w_p1 = kwargs.get('w_p1', 0)
        self.w_p2 = kwargs.get('w_p2', 0)
        self.fcn = 0
        self.GanomT = 1
        self.tf2od = 0
        self.G_bp = 0
        self.tipo = tipo
        self.w = 0
        self.amp = 0
        self.fase = 0
        self.wc = 0.0
        self.Rl = 0
        self.top = ''
        self.roots = 0
        self.den_norm = 0
        self.Bw = 0
        self.Bp = 0
        self.Bs = 0
        self.criterio = 0

    def ganho_bp(self, G):
        if(G >= self.a_p):
            self.G_bp = -G
            self.a_s = self.a_s + self.G_bp
            self.a_p = self.a_p + self.G_bp
        else:
            raise ValueError("A atenuação na banda de passagem deve ser menor que o ganho na banda de passagem")

    def ordem(self, **kwargs):
        result = 0
        if self.tipo == 'hp':
            result = int(np.ceil(np.log10((np.power(
            ((np.power(10.0, (-self.a_s / 10.0)) - 1.0) / (np.power(10.0, (-self.a_p / 10.0)) - 1.0)),
            0.5))) / np.log10((self.w_p / self.w_s))))
        if self.tipo == 'lp':
            result = int(np.ceil(np.log10((np.power(
            ((np.power(10.0, (-self.a_s / 10.0)) - 1.0) / (np.power(10.0, (-self.a_p / 10.0)) - 1.0)),
            0.5))) / np.log10((self.w_s / self.w_p))))
        if(self.tipo == 'bp'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.log10((np.power(
            ((np.power(10.0, (-self.a_s / 10.0)) - 1.0) / (np.power(10.0, (-self.a_p / 10.0)) - 1.0)),
            0.5))) / np.log10((self.Bs / self.Bp))))
        if(self.tipo == 'bs'):
            self.Bp = self.w_p2 - self.w_p1
            self.Bs = self.w_s2 - self.w_s1
            result = int(np.ceil(np.log10((np.power(
            ((np.power(10.0, (-self.a_s / 10.0)) - 1.0) / (np.power(10.0, (-self.a_p / 10.0)) - 1.0)),
            0.5))) / np.log10((self.Bp / self.Bs))))
        result = kwargs.get('ord', result)
        return result

    def fq_corte(self, **kwargs):
        result = None
        if ('ordem' in kwargs):
            ordem = kwargs['ordem']
        else:
            ordem = self.ordem()
        if(self.tipo == 'hp'):
            result = self.w_p*(np.power((np.power(10.0, ((-1.0) * self.a_p / 10.0)) - 1.0), (1.0 / (2.0*ordem))))
        elif(self.tipo == 'lp'):
            result = self.w_p / (np.power((np.power(10.0, ((-1.0) * self.a_p / 10.0)) - 1.0), (1.0 / (2.0 * ordem))))
        elif(self.tipo == 'bp'):
            result = np.sqrt(self.w_p1*self.w_p2)
        elif(self.tipo == 'bs'):
            result = np.sqrt(self.w_s1*self.w_s2)
        self.wc = result
        return result

    def raizes_normal(self, **kwargs):
        ordem = kwargs.get('ordem', self.ordem())
        S_k = np.zeros(ordem, dtype=complex)
        for i in range(1, ordem + 1):
            S_k[i - 1] = -np.sin(np.pi * (2 * i - 1) / (2 * ordem)) + (
                        1j * np.cos(np.pi * (2 * i - 1) / (2 * ordem)))  # cálculo das raízes
        self.roots = S_k
        return self.roots

    def transfunc(self, polos, **kwargs):
        wc = kwargs.get('wc', 0)
        w0 = kwargs.get('w0', 0)
        Bw = kwargs.get('bw', self.Bp)
        ordem = kwargs.get('ord', self.ordem)
        resp = kwargs.get('response', self.tipo)
        G_db = kwargs.get('G', self.G_bp)
        fcn = 0
        if resp == 'lp':
            self.den_norm = np.real(np.poly(polos))
            denm = np.zeros(len(self.den_norm))
            for i in range(0, len(polos) + 1):
                denm[i] = self.den_norm[i] * np.power(wc, i)
            num = denm[-1]
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*num, denm)
        if resp == 'hp':
            self.den_norm = np.real(np.poly(polos))
            denm = np.zeros(len(polos) + 1)
            for i in range(0, len(polos) + 1):
                denm[i] = self.den_norm[len(polos) - i] * np.power(wc, i)
            num = np.zeros(len(polos) + 1)
            num[0] = denm[0]
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*num, denm)
        if(resp == 'bp'):
            self.den_norm = np.real(np.poly(polos))
            [num, den] = signal.lp2bp(self.den_norm[-1], self.den_norm, w0, Bw)
            print('EMBAIXO')
            print(self.den_norm)
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*num, den)
        if(resp == 'bs'):
            self.den_norm = np.real(np.poly(polos))
            [num, den] = signal.lp2bs(self.den_norm[-1], self.den_norm, w0, Bw)
            fcn = signal.TransferFunction((pow(10, -G_db/20.0))*num, den)
        self.fcn = fcn
        return fcn

    def transfunc2(self, t):
        ordem = self.ordem()
        polos2o = np.zeros(2)
        tf2o = []
        if(t.lower() == 'sk' and (self.tipo == 'lp' or self.tipo == 'hp' )):
            if(ordem%2 == 0):
                for k in range(0, int(ordem/2)):
                    polos2o = [self.roots[k], self.roots[-k-1]]
                    tf = self.transfunc(polos2o, wc = self.wc, G = 0)
                    tf2o.append(tf)
            else:
                for k in range(0, int(ordem/2)+1):
                    if(k != int(ordem/2)):
                        polos2o = [self.roots[k], self.roots[-k-1]]
                        tf = self.transfunc(polos2o, wc = self.wc, G = 0)
                        tf2o.append(tf)
                    else:
                        polo1o = [self.roots[int(ordem/2)]]
                        tf = self.transfunc(polo1o, wc = self.wc, G = self.G_bp)
                        tf2o.append(tf)
        if(t.lower() == 'notch' and (self.tipo == 'bs')):
            polosNo = np.roots(self.fcn.den)
            zerosNo = np.roots(self.fcn.num)
            modZ = np.abs(zerosNo)
            modP = np.abs(polosNo)
            polosNo = polosNo[np.argsort(modP)]
            zerosNo = zerosNo[np.argsort(modZ)]
            for k in range(0, int(ordem)):
                zeros2o = [zerosNo[2*k], zerosNo[2*k + 1]]
                polos2o = [polosNo[2*k], polosNo[2*k + 1]]
                print(self.fcn.den)
                print(self.fcn.num)
                num = np.poly((1j)*np.imag(zeros2o))
                den = np.poly(polos2o)
                tf = signal.TransferFunction(num, den)
                tf2o.append(tf)
        if(t.lower() == 'friend' and self.tipo == 'bp'):
            print(self.fcn.den)
            polosNo = np.roots(self.fcn.den)
            tf2o = []
            Qv = []
            num = [1, 0]
            for i in range(ordem):
                polos2o = [polosNo[2*i], np.conjugate(polosNo[2*i])]
                den = np.poly(polos2o)
                tf_aux = signal.TransferFunction(num, den)
                tf2o.append(tf_aux)
                Q_aux = np.sqrt(tf_aux.den[-1])/tf_aux.den[-2]
                Qv.append(Q_aux)
        self.tf2od = tf2o
        return tf2o, Qv


    def elementosAtivos(self, t, Q, **kwargs):
        Comp = {}
        w0 = kwargs.get('w0', 1)
        if(t.lower() == 'friend'):
            Cn = (1e-6)*float(input('Valor do capacitor fixo (uF): '))
            Ki = 1/(2*Q*w0*Cn) 
            Comp["C1"] = Cn
            Comp["C2"] = Cn
            Comp["R1"] = 1*Ki
            Comp["R2"] = (4*pow(Q, 2))*Ki
        return Comp

    def graphpoints(self, fcn, max_f, min_f, points):
        dt = float((max_f - min_f)/points)
        self.w, self.amp, self.fase = fcn.bode(w=np.arange(min_f, max_f, dt))
        if(self.tipo == 'hp'  or self.tipo == 'lp'):
            kp = np.int(np.round((self.w_p - min_f)/dt))
            kr = np.int(np.round((self.w_s - min_f)/dt))
            if(self.amp[kp] >= (self.a_p - self.G_bp) and self.amp[kr] <= (self.a_s - self.G_bp)):
                print('Pontos de projeto ATENDIDOS!')
                self.criterio = 0               #Não há como otimizar, apenas aumentar a ordem
            else:
                print('Pontos de projeto NÃO ATENDIDOS!')
                self.criterio = 0               #Não há como otimizar, apenas aumentar a ordem
        if(self.tipo == 'bp' or self.tipo == 'bs'):
            kp1 = np.int(np.round((self.w_p1 - min_f)/dt))-1
            kp2 = np.int(np.round((self.w_p2 - min_f)/dt))-1
            kr1 = np.int(np.round((self.w_s1 - min_f)/dt))-1
            kr2 = np.int(np.round((self.w_s2 - min_f)/dt))-1
            if(self.amp[kp1] >= (self.a_p - self.G_bp)) and (self.amp[kp2] >= (self.a_p - self.G_bp)) and (self.amp[kr1] <= (self.a_s - self.G_bp)) and (self.amp[kr2] <= (self.a_s - self.G_bp)):
                print('Pontos de projeto: ATENDIDOS!')
                self.criterio = 0
            else:
                print('Pontos de projeto: NÃO ATENDIDOS!')
                self.criterio = -1


    def plot_bode(self, fcn, **kwargs):
     
        fig, ax = plt.subplots()  # cria os plots
        ax.semilogx(self.w, self.amp)  # gráfico do tipo semilog
        ax.set(xlabel="Frequência (rad/s)", ylabel="Amplitude em dB",
               title="Resposta em amplitude (BTRWRTH n = %d)" %self.ordem()) # configuração de plot label
        if(self.tipo == 'lp' or self.tipo == 'hp'):
            bp = ax.scatter(self.w_p, self.a_p - self.G_bp)  # ponto de projeto de passagem
            br = ax.scatter(self.w_s, self.a_s - self.G_bp)  # ponto de projeto de rejeição
            ax.legend((bp, br), ("P. Projeto (passagem)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        elif(self.tipo == 'bp' or self.tipo == 'bs'):
            bp1 = ax.scatter(self.w_p1, self.a_p - self.G_bp)
            bp2 = ax.scatter(self.w_p2, self.a_p - self.G_bp)  # ponto de projeto de passagem
            bs1 = ax.scatter(self.w_s1, self.a_s - self.G_bp)  # ponto de projeto de rejeição
            bs2 = ax.scatter(self.w_s2, self.a_s - self.G_bp)
            ax.legend((bp1, bp2, bs1, bs2), ("P. Projeto (passagem)", "P. Projeto (passagem)", "P. Projeto (rejeição)", "P. Projeto (rejeição)"), loc='lower left', fontsize=8)
        ax.margins(x=0)
        ax.margins(y=0.05)  # margem y
        ax.grid(True, which="both")  # grid
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        return fig

    def elements(self, wc, topologia, R,**kwargs):                  #Calcula os elementos escalonados em frequência
        if ('ordem' in kwargs):                                     #Calcula a ordem, caso seja chamado antes
            ordem = kwargs['ordem']
        else:
            ordem = self.ordem()
        bw = kwargs.get('bw', 0)
        self.top = topologia
        self.Rl = R
        nomes = np.empty(ordem, dtype="<U2")
        val = np.zeros(ordem)
        elem_table = 2*abs(np.real(self.roots))                     #Utiliza a formula Cauer para calcular os elementos da tabela Butterworth
        esq_dir = elem_table[::-1]                                  #Coloca os elementos na ordem inversa da tabela, para coincidir o primeiro elemento a depender da self.topologia
        if(self.tipo == 'lp'):                                      #Cálculo para FPB
            elem_fpb = esq_dir/wc                                   #Escalonamento em frequência (independe se é L ou C)
            if(topologia == 'clc'):                                 #Inicia com capacitor
                val, nomes = self.__Zscale(ordem, elem_fpb, R, topologia)
            if(topologia == 'lcl'):
                val, nomes = self.__Zscale(ordem, elem_fpb, R, topologia)
        if(self.tipo == 'hp'):
            elem_fpa = (np.float_power(esq_dir, -1))*(1/wc)                                   #Escalonamento em frequência (independe se é L ou C)
            if(topologia == 'clc'):
                val, nomes = self.__Zscale(ordem, elem_fpa, R, topologia)
            if(topologia == 'lcl'):
                val, nomes = self.__Zscale(ordem, elem_fpa, R, topologia)
        if(self.tipo == 'bp' or self.tipo == 'bs'):
            val, nomes = self.__Zscale(ordem, esq_dir, R, topologia, bw = bw)
        return val, nomes                                           #Retorna os elementos

    def __Zscale(self, N, elementos, R, topologia, **kwargs):
        bw = kwargs.get('bw', 0)
        name = np.empty(2*N, dtype="<U2")
        if(self.tipo == 'hp' or self.tipo == 'lp'):
            elem_final2 = np.zeros(N)
        elif(self.tipo == 'bs' or self.tipo == 'bp'):
            elem_final2 = np.zeros(2*N)
        if(topologia == 'clc'):                                     #Primeiro elemento é capacitor
            for i in range(1, N+1):                                 #Varre todos os elementos
                if(i%2 != 0):                                       #Se começa com capacitor, todos os elementos ímpares são capacitores
                    elem_final2[i-1] = elementos[i-1]/R              #Capacitor
                    name[i-1] = "C" + str(N-i+1)                    #Nome do elemento + numero (ordem do número segue o Paarman L.)
                else:                                               #Se par, é indutor
                    elem_final2[i-1] = elementos[i-1]*R              #Indutor
                    name[i-1] = "L" + str(N-i+1)                    #Nome do elemento + numero (ordem do número segue o Paarman L.)
        if(topologia == 'lcl'):                                     #Primeiro elemento é indutor
            for i in range(1, N+1):                                 #Varre todos os elementos
                if(i%2 != 0):                                       #Se começa com indutor, todos os elementos ímpares são indutores
                    elem_final2[i-1] = elementos[i-1]*R              #Indutor
                    name[i-1] = "L" + str(N-i+1)                    #Nome do elemento + numero (ordem do número segue o Paarman L.)

                else:                                               #Se par, então é capacitor
                    elem_final2[i-1] = elementos[i-1]/R              #Capacitor
                    name[i-1] = "C" + str(N-i+1)                    #Nome do elemento + numero (ordem do número segue o Paarman L.)
        if(topologia == 's' and self.tipo == 'bp'):
            for i in range(1, N+1):
                if(i%2 != 0):
                    elem_final2[2*i-2] = elementos[i-1]*(R/bw)
                    elem_final2[2*i-1] = bw*(1/(elementos[i-1]*np.power(self.wc,2)*R))
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
                else:
                    elem_final2[2*i-2] = elementos[i-1]/(bw*R)
                    elem_final2[2*i-1] = (bw*R)*(1/(elementos[i-1]*np.power(self.wc,2)))
                    name[2*i-2] = 'C' + str(N-i+1)
                    name[2*i-1] = 'L' + str(N-i+1)
        if(topologia == 'p' and self.tipo == 'bp'):
            for i in range(1, N+1):
                if(i%2 == 0):
                    elem_final2[2*i-2] = elementos[i-1]*(R/bw)
                    elem_final2[2*i-1] = bw*(1/(elementos[i-1]*np.power(self.wc,2)*R))
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
                else:
                    elem_final2[2*i-2] = elementos[i-1]/(bw*R)
                    elem_final2[2*i-1] = (bw*R)*(1/(elementos[i-1]*np.power(self.wc,2)))
                    name[2*i-2] = 'C' + str(N-i+1)
                    name[2*i-1] = 'L' + str(N-i+1)
        if(topologia == 's' and self.tipo == 'bs'):
            for i in range(1, N+1):
                if(i%2 != 0):
                    elem_final2[2*i-2] = R/(bw*elementos[i-1])
                    elem_final2[2*i-1] = (bw*(elementos[i-1]))/(R*np.power(self.wc,2))
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
                else:
                    elem_final2[2*i-2] = R*(bw*elementos[i-1])/(np.power(self.wc,2))
                    elem_final2[2*i-1] = 1/(bw*elementos[i-1]*R)
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
        if(topologia == 'p' and self.tipo == 'bs'):
            for i in range(1, N+1):
                if(i%2 == 0):
                    elem_final2[2*i-2] = R/(bw*elementos[i-1])
                    elem_final2[2*i-1] = (bw*(elementos[i-1]))/(R*np.power(self.wc,2))
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
                else:
                    elem_final2[2*i-2] = R*(bw*elementos[i-1])/(np.power(self.wc,2))
                    elem_final2[2*i-1] = 1/(bw*elementos[i-1]*R)
                    name[2*i-2] = 'L' + str(N-i+1)
                    name[2*i-1] = 'C' + str(N-i+1)
        return elem_final2, name

    def exporttxt(self, arq_nome):
        with open(arq_nome, 'w') as f:
            f.write('Detalhamento do filtro\n')
            f.write('Pontos de projeto escolhidos:\n')
            f.write('   Filtro: Butterworth\n')
            f.write('   Resposta: ' + str(self.tipo.upper()) + '\n')
            f.write('   Ordem: ' + str(self.ordem()) + '\n')
            f.write('   Frequência de corte: ' + str(self.wc) + '\n')
            f.write('\nComponentes do circuito:\n')
            f.write('   Topologia: ' + str(self.top) + '\n')
            valor, nome = self.elements(self.wc, self.top, self.Rl)
            f.write("   Elementos escalonados: \n")
            for i in range(0, self.ordem()):
                f.write('     %s = %s\n' %(nome[i], eng_string(valor[i], si=True)))
            f.write("\n")         
       
    def testa_par(self):
        if (self.tipo == 'hp') and ((self.w_s >= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um HP")
        if (self.tipo == 'lp') and ((self.w_s <= self.w_p) or (self.a_s >= self.a_p)):
            raise ValueError("Inversão dos parâmetros para um LP")

    def bw_increase(self, pct):
        if(self.tipo == 'bs'):
            pct = pct*(-1)
        Bw_mod = (self.w_p2-self.w_p1)*(1+pct/100)
        tf1 = self.transfunc(self.roots, w0 = self.wc, bw = Bw_mod)
        return tf1, Bw_mod

def eng_string( x, format='%s', si=False):
    sign = ''
    if x < 0:
        x = -x
        sign = '-'
    exp = int( math.floor( math.log10( x)))
    exp3 = exp - ( exp % 3)
    x3 = x / ( 10 ** exp3)

    if si and exp3 >= -24 and exp3 <= 24 and exp3 != 0:
        exp3_text = 'yzafpnum kMGTPEZY'[ int(( exp3 - (-24)) / 3)]
    elif exp3 == 0:
        exp3_text = ''
    else:
        exp3_text = 'e%s' % exp3

    return ( '%s'+format+'%s') % ( sign, "{:.3f}".format(x3), exp3_text)
