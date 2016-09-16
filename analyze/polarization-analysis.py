import csv
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

def Read_File(filename):
	tl = []
	ptl = []
	with open(filename) as csvfile:
		reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
		counter = 0
		for row in reader:
			counter += 1
			x = float(row[0])
			y = float(row[1])
			tl.append(x)
			ptl.append(y)

	n = len(tl) / 2
	for i in xrange(n):
		tl.pop(0)
		ptl.pop(0)

	t = np.asarray(tl)
	pt = np.asarray(ptl)
	return t, pt


def Auto_Corr(t,pt):
	n = t.size
	m = n / 2
	di = m / 20
	ave = pt.mean()
	var = pt.var()
	taui = 0
	ci = 0
	taul = []
	cl = []
	for i in xrange(0,m,di):
		ci = 0
		for j in xrange(n-i):
			ci += (pt[j] - ave)*(pt[j+i] - ave) / (n - i)
		taui = t[i] - t[0]
		ci /= var
		taul.append(taui)
		cl.append(ci)
	tau = np.asarray(taul)
	c = np.asarray(cl)
	ic = 0
	while (c[ic] > 0.01):
		ic += 1
	ic *= di
	return tau,c,ic


def Find_Polarity(t,pt,ic):
	n = pt.size / ic
	error = math.sqrt(pt.var() / n)
	ave = pt.mean()
	return ave,error
	



arg_list = sys.argv[1:]
for i in arg_list:
	filename = i
	t,pt = Read_File(filename)
	start = filename.find('Dr=') + 2
	end = filename.find('-K=', start)
	tau,c,ic = Auto_Corr(t,pt)
#	for j in xrange(tau.size):
#		print tau[j],c[j]
	ave, error = Find_Polarity(t,pt,ic)
	print filename[start+1:end],ave,error




