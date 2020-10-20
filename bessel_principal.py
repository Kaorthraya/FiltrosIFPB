from bessel import f_bessel
import numpy as np
import matplotlib.pyplot as plt

resposta = 'hp'
wp = (26000.0)*(2.0*np.pi)         #Frequência de passagem
ws = (4000.0)*(2.0*np.pi)          #Frequência de rejeição
ap = 20*np.log10(0.9)              #Valor de amplitude na banda de passagem
a_s = 20*np.log10(0.1)             #Valor de amplitude na banda de rejeição
ordem = 1                          #Ordem inicial
ordem_max = 20                     #Ordem máxima do filtro (até obter convergência ou falhar)

filtro = f_bessel(resposta, a_p = ap, a_s = a_s, w_p = wp, w_s = ws)

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
            kf = wp/wpl
        elif(resposta == 'hp'):
            kf = wp*wpl        
        tf_transf1 = filtro.transf(ordem, resposta, kf = kf)
        filtro.points(1e3, 1e8, 1e5)            #Gera o gráfico escalonado na frequẽncia correta
        filtro.test_point()                     #Testa o ponto de rejeição
        ordem = ordem+1                         #Itera a ordem

if(filtro.criterio == 0):
    print('\nCONVERGIU!')
    print('Ordem aceita: %d' %(ordem-1))

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
    if(resposta == 'hp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(wp, filtro.w[-1], alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(0, ws, alpha=0.1, color='red')
    elif(resposta == 'lp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(0, wp, alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(ws, filtro.w[-1], alpha=0.1, color='red')
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
    if(resposta == 'hp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(wp, filtro.w[-1], alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(0, ws, alpha=0.1, color='red')
    elif(resposta == 'lp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(0, wp, alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(ws, filtro.w[-1], alpha=0.1, color='red')
    #plt.savefig('fase.png', dpi=600)

    #Configuração de plotagem do Group delay
    fig3 = plt.figure(3)
    ax = fig3.gca()
    filtro.gp_delay()
    ax.semilogx(filtro.w, filtro.gpdelay)
    ax.margins(x=0)
    ax.margins(y=0.05)
    ax.set_title('Delay de grupo')
    ax.set_ylabel('Delay em segundos')
    ax.set_xlabel('Frequência em rad/s')
    ax.grid(True, which="both")
    if(resposta == 'hp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(wp, filtro.w[-1], alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(0, ws, alpha=0.1, color='red')
    elif(resposta == 'lp'):
        plt.axvline(wp, color = 'green')
        ax.axvspan(0, wp, alpha=0.1, color='green')
        plt.axvline(ws, color = 'red')
        ax.axvspan(ws, filtro.w[-1], alpha=0.1, color='red')
    #plt.savefig('gd.png', dpi=600)
    plt.show()
else:
    print('FIltro NÃO convergiu')
