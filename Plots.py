import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolors
import scipy 
from scipy.io.wavfile import write 

#### CUERDA 

## EXTREMOS FIJOS 

x = np.loadtxt('x_x.dat')
phi_0 = np.loadtxt('phi_0.dat')
phi_12 = np.loadtxt('phi_12.dat')
phi_25 = np.loadtxt('phi_25.dat')
phi_50 = np.loadtxt('phi_50.dat')

phi = [phi_0, phi_12, phi_25, phi_50]

for phi_i in phi:
	plt.plot(x, phi_i)

plt.savefig("cuerda_extremos_fijos.png")
#plt.show()
plt.clf()

## PERTURBACION EN UN EXTREMO

psi_0 = np.loadtxt('psi_0.dat')
psi_12 = np.loadtxt('psi_12.dat')
psi_25 = np.loadtxt('psi_25.dat')
psi_50 = np.loadtxt('psi_50.dat')

psi = [psi_0, psi_12, psi_25, psi_50]

for psi_i in psi:
	plt.plot(x, psi_i)

plt.xlim(0.0, 0.2) 

plt.savefig("cuerda_perturbacion_periodica.png")
#plt.show()
plt.clf()

#### ARCHIVO DE SONIDO 

phi_med = np.loadtxt('phi_med.dat')

scipy.io.wavfile.write('sonido.WAV', phi_med.shape[0], phi_med)
#print(scipy.io.wavfile.read("sonido.wav"))

#### TAMBOR 

xi_x0 = np.loadtxt('xi_x0.dat')
xi_y0 = np.loadtxt('xi_y0.dat')
xi_x12 = np.loadtxt('xi_x12.dat')
xi_y12 = np.loadtxt('xi_y12.dat')
xi_x25 = np.loadtxt('xi_x25.dat')
xi_y25 = np.loadtxt('xi_y25.dat')
xi_x50 = np.loadtxt('xi_x50.dat')
xi_y50 = np.loadtxt('xi_y50.dat')

t_0 = np.vstack([xi_x0 ,xi_y0])
t_12 = np.vstack([xi_x12 ,xi_y12])
t_25 = np.vstack([xi_x25 ,xi_y25])
t_50 = np.vstack([xi_x50 ,xi_y50])

t = [t_0, t_12, t_25, t_50]

tambor = plt.figure()

Im_t_0 = tambor.add_subplot(221)
Im_t_12 = tambor.add_subplot(222)
Im_t_25 = tambor.add_subplot(223)
Im_t_50 = tambor.add_subplot(224)

Im_t = [Im_t_0, Im_t_12, Im_t_25, Im_t_50]

for Im_t_i, t_i in zip(Im_t, t):
	Im_t_i.imshow(t_i)
	Im_t_i.axis('off')
	Im_t_i.set_ylim(0, 100)

plt.savefig("membrana_tambor.png")
#plt.show()
plt.clf()


