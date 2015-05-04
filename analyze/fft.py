import numpy
import math

a = numpy.fromfile("pair-dist.dat", float,-1," ")

N = int(math.sqrt(len(a)/3))

v = numpy.zeros((N,N))

L = -a[0]

for i in xrange(N):
	for j in xrange(N):
		v[i,j] = a[i*3*N+3*j+2]

fv = numpy.fft.fft2(v) * ((2*L/N)**2)
fvnorm = numpy.absolute(fv)

f = file("pair-dist-fft.dat",'w')

for i in xrange(N):
	for j in xrange(N):
		a[i*3*N+3*j+2] = fvnorm[i,j]
		if (j < N/2):
			a[i*3*N+3*j+1] = (math.pi*j / L) / (2*math.pi)
		else:
			a[i*3*N+3*j+1] = (math.pi*(j-N) / L) / (2*math.pi)
		if (i < N/2):
			a[i*3*N+3*j] = (math.pi*i / L) / (2*math.pi)
		else:
			a[i*3*N+3*j] = (math.pi*(i-N) / L) / (2*math.pi)

for i in xrange(N/2,N):
	for j in xrange(N/2,N):
		f.write(str(a[i*3*N+3*j]) + "\t" + str(a[i*3*N+3*j+1]) + "\t" + str(a[i*3*N+3*j+2]) + "\n")
	for j in xrange(N/2):
		f.write(str(a[i*3*N+3*j]) + "\t" + str(a[i*3*N+3*j+1]) + "\t" + str(a[i*3*N+3*j+2]) + "\n")
	f.write("\n")

for i in xrange(N/2):
	for j in xrange(N/2,N):
		f.write(str(a[i*3*N+3*j]) + "\t" + str(a[i*3*N+3*j+1]) + "\t" + str(a[i*3*N+3*j+2]) + "\n")
	for j in xrange(N/2):
		f.write(str(a[i*3*N+3*j]) + "\t" + str(a[i*3*N+3*j+1]) + "\t" + str(a[i*3*N+3*j+2]) + "\n")
	f.write("\n")



