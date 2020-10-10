from bessel import *
import matplotlib.pyplot as plt

resposta = 'hp'
wp = 5000.0         #Frequência de passagem
ws = 1000.0         #Frequência de rejeição
ap = -1.0           #Valor de amplitude na banda de passagem
a_s = -20.0         #Valor de amplitude na banda de rejeição
ordem = 1           #Ordem inicial
ordem_max = 10      #Ordem máxima do filtro (até obter convergência ou falhar)

filtro = f_bessel(resposta, a_p = ap, a_s = a_s, w_p = wp, w_s = ws)

while(ordem <= ordem_max and filtro.criterio == -1):
        print('\nOrdem: %d' %ordem)
        den = filtro.coeff(ordem)
        tf = filtro.tf_norm()
        filtro.points(0.01, 5, 10e3)            #Gera pontos pro gráfico não normalizado
        ind, err = filtro.wp_norm(0.001, 100)   #seta ordem e tolerância do Pégaso
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
        filtro.points(1,10e4,100e3)
        filtro.test_point()
        ordem = ordem+1

if(filtro.criterio == 0):
    print('CONVERGIU!')
    print('Ordem aceita: %d' %(ordem-1))

    fig1 = plt.figure(1)
    ax = fig1.gca()
    ax.semilogx(filtro.w, filtro.amp)
    ax.margins(x=0)
    ax.margins(y=0.05)  
    ax.grid(True, which="both")
    ax.set_title('Magnitude')
    ax.set_ylabel('Magnitude em dB')
    ax.set_xlabel('Frequência em rad/s')
    filtro.plot_scatter(ax)

    fig2 = plt.figure(2)
    ax = fig2.gca()
    ax.semilogx(filtro.w, filtro.fase)
    ax.margins(x=0)
    ax.margins(y=0.05) 
    ax.set_title('Fase')
    ax.set_ylabel('Fase em graus')
    ax.set_xlabel('Frequência em rad/s')
    ax.grid(True, which="both")

    plt.show()
else:
    print('FIltro NÃO convergiu')
