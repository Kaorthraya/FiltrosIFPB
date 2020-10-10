import numpy as np
from scipy import signal

class f_bessel:
    def __init__(self, tipo, **kwargs):
        self.a_p = kwargs['a_p']
        self.a_s = kwargs['a_s']
        self.w_s = kwargs.get('w_s', 0)
        self.w_p = kwargs.get('w_p', 0)
        self.w_s1 = kwargs.get('w_s1', 0)
        self.w_s2 = kwargs.get('w_s2', 0)
        self.w_p1 = kwargs.get('w_p1', 0)
        self.w_p2 = kwargs.get('w_p2', 0)
        self.criterio = -1
        self.min_f = 0
        self.max_f = 0
        self.step = 0
        self.tf = 0
        self.coeficientes = 0
        self.w = 0
        self.amp = 0
        self.fase = 0
        self.tipo = tipo
        self.gpdelay = 0

    def coeff(self, N):
        aux = np.zeros(N+1, dtype=float)
        aux[0] = 1.0
        for i in range(0, N):
            aux[i+1] = (2*(N-i)*aux[i])/((2*N-i)*(i+1))
        self.coeficientes = aux[::-1]
        return self.coeficientes

    def tf_norm(self):
        self.tf = signal.TransferFunction(1, self.coeficientes)
        return self.tf

    def points(self, min_f, max_f, num_points):
        self.min_f = min_f
        self.max_f = max_f
        self.step = float((max_f - min_f)/(num_points))
        self.w, self.amp, self.fase = self.tf.bode(w = np.arange(min_f, max_f, self.step))
        self.fase = self.fase*(np.pi/180.0)

    def wp_norm(self, tol, itmax):
        #Algoritmo Pégaso
        amp_aux = self.amp - self.a_p
        a = self.min_f+0.01
        b = self.max_f-0.01
        Fa = amp_aux[self.index(a)]
        Fb = amp_aux[self.index(b)]
        i = 0
        while(i <= itmax):
            dx = (-Fb*(b-a))/(Fb - Fa)
            x = b+dx
            Fx = amp_aux[self.index(x)]
            i = i+1
            if(abs(dx) <= tol and abs(Fx) <= tol):
                return self.index(x), False
            if(Fx * Fb < 0):
                a = b
                Fa = Fb
            else:
                Fa = (Fa* Fb)/(Fb+Fx)
            b = x
            Fb = Fx
        return x, True

    def transf(self, ordem, tipo, **kwargs):
        kf = kwargs.get('kf', 0)
        Bw = kwargs.get('bw', 0)
        if(tipo == 'lp'):
            den = np.zeros(ordem+1, dtype=float)
            for i in range(0, ordem+1):
                den[i] = self.coeficientes[i]*pow(kf, i)
            num = den[-1]
            self.tf = signal.TransferFunction(num, den)
        if(tipo == 'hp'):
            num = np.zeros(ordem+1, dtype=float)
            den = np.zeros(ordem+1, dtype=float)
            num, den = signal.lp2hp(1, self.coeficientes, kf)
            self.tf = signal.TransferFunction(num, den)
        return self.tf
    
    def test_point(self):
        if(self.tipo == 'lp' or self.tipo == 'hp'):
            if(self.amp[self.index(self.w_s)] <= self.a_s):
                print('Ponto de rejeição aceito!')
                self.criterio = 0
            else:
                print('Ponto de rejeição NÃO aceito!')
                self.criterio = -1

    def plot_scatter(self, ax):
        if(self.tipo == 'lp' or self.tipo == 'hp'):
            bp = ax.scatter(self.w_p, self.a_p, color = 'green')  # ponto de projeto de passagem
            br = ax.scatter(self.w_s, self.a_s, color = 'red')  # ponto de projeto de rejeição
            ax.legend((bp, br), ("P. Projeto (passagem)", "P. Projeto (rejeição)"), loc='lower right', fontsize=8)

    def gp_delay(self):
        gd = (-1)*(np.diff((self.fase))/np.diff(self.w))
        self.gpdelay = np.append(gd, 0)  

    def index(self, freq):
        ind = np.int(np.round((freq - self.min_f)/self.step))
        return ind
