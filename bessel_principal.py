from bessel import *
import numpy as np
from printft import printTF
import matplotlib.pyplot as plt

print('Cálculo de filtros Bessel\n')
print('Opções de resposta:\n[lp] - lowpass\n[hp] - highpass\n[bp] - bandpass\n[bs] - bandstop\n')
resposta = input('Digite a resposta do filtro: ')
if(resposta.lower() == 'lp' or resposta.lower() == 'hp'):
    ap_dB = float(input('Valor da atenuação na banda de passagem: '))
    as_dB = float(input('Valor da atenuação na banda de rejeição: '))
    Wp = float(input('Valor de frequência (passagem): '))
    Ws = float(input('Valor de frequência (rejeição): '))
    filtro = f_bessel(resposta.lower(),a_p = ap_dB,
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
    filtro = f_bessel(resposta.lower(),a_p = ap_dB,
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

ordem = 1                          #Ordem inicial
ordem_max = 20                     #Ordem máxima do filtro (até obter convergência ou falhar)

while(ordem <= ordem_max and filtro.criterio == -1):
        print('\nOrdem: %d' %ordem)
        den = filtro.coeff(ordem)
        tf = filtro.tf_norm()                   #Gera os coeficientes de Bessel
        filtro.points(0.01, 5, 10e3)            #Gera pontos pro gráfico não normalizado
        ind, err = filtro.wp_norm(0.01, 100)   #seta ordem e tolerância do Pégaso
        if(err == True):
            break
        else:
            wpl = filtro.w[ind]
            print('Pégaso convergiu, Wp encontrado para ordem: %d' %ordem)
            print('Frequência Wp: %2.2f' %wpl)
        if(resposta == 'lp'):
            kf = Wp/wpl
            tf_transf1 = filtro.transf(ordem, resposta, G = Ganho,  kf = kf)
        elif(resposta == 'hp'):
            kf = Wp*wpl
            tf_transf1 = filtro.transf(ordem, resposta, G = Ganho, kf = kf)
        elif(resposta == 'bp'):
            w0 = np.sqrt(Wp1*Wp2)
            Bw = Wp2 - Wp1
            tf_normalizada = filtro.transf(ordem, 'lp', kf = 1/wpl)
            tf_transf1 = filtro.transf(ordem, resposta, w0 = w0, bw = Bw, G = Ganho)
        printTF(tf_transf1.num, tf_transf1.den)
        filtro.points(10, 1e6, 1e6)  
        filtro.test_point()                     #Testa o ponto de rejeição
        if(filtro.criterio == 0):
            break
        else:
            ordem = ordem+1

if(filtro.criterio == 0):
    print('\nCONVERGIU!')
    print('Ordem aceita: %d' %(ordem))

    #Configurações do SOS
    Tfs, Q = filtro.transfunc2(2*ordem, 'friend')
    for i in range(0, len(Tfs)):
        printTF(Tfs[i].num, Tfs[i].den)
        print('Valor de Q: %2.3f' %Q[i])
        componentes = filtro.elementosAtivos('friend', Q[i], w0 = np.sqrt(Tfs[i].den[-1]))
        for x in componentes:
            print('%s: %s' %(x, eng_string(componentes[x], '%s', si = True)))

    #Configurações de plotagem da magnitude
    fig1 = plt.figure(1)
    ax = fig1.gca()
    ax.semilogx(filtro.w, filtro.amp)
    ax.margins(x=0)
    ax.margins(y=0.05)  
    ax.grid(True, which="both")
    ax.set_title('Magnitude')
    ax.set_ylabel('Magnitude em dB')
    ax.set_xlabel('Frequência em rad/s') 
    ax.set_ylim(-100, 5)
    filtro.plot_scatter(ax)
    #plt.savefig('mag.png', dpi=600)
    
    #Configuração de plotagem da fase
    fig2 = plt.figure(2)
    ax = fig2.gca()
    ax.semilogx(filtro.w, filtro.fase)
    ax.margins(x=0)
    ax.margins(y=0.05) 
    ax.set_title('Fase')
    ax.set_ylabel('Fase em radianos')
    ax.set_xlabel('Frequência em rad/s')
    ax.grid(True, which="both")
    #plt.savefig('fase.png', dpi=600)
    plt.show()
else:
    print('FIltro NÃO convergiu')
