import numpy as N
import copy
import sys
#import pylab 

def three_magn_corr(configs,vec):#<psi_ground|s_i^+ s_j^+ s_o^+ |psi_3magGround>
	L=len(configs[0])
	corr=N.zeros((L,L),dtype=float)
	special=(L/2)-1
	for h, config in enumerate(configs):
		for i in range(L):
			for j in range(L):
				if (config[special]==1):
					corr[i,j]+=0
				if (i!=special) and (j!=special) and (i!=j):
					if (config[special]==0) and (config[i]==config[j]==0):
						corr[i,j]+=8*vec[h]*vec[h]
				if (i==j):
					if (config[special]==0) and (config[i]==-1):
						corr[i,j]+=8*vec[h]*vec[h]
				if (j==special):
					if (config[special]==-1) and (config[i]==0):
						corr[i,j]+=8*vec[h]*vec[h]
				if (i==special):
					if (config[special]==-1) and (config[j]==0):
						corr[i,j]+=8*vec[h]*vec[h]
	return corr

def search_config(config,configs):
	present= False
	L=len(configs[0])
	for loc,c in enumerate(configs):
		match=True
		for n in range(L):
			if (configs[loc][n]!=config[n]):
				match=False
				n=L-1
		if (match==True):
			return loc,True
	return -1,False

L=int(sys.argv[1])
'''
Bkg=float(sys.argv[2])
g=2.14
JFMK=14.8
DK=+5.2
DK=2.0*DK
print "Z DIAG "
print "L        = ",L
print "JFM (K)  = ",JFMK 
print "D   (K)  = ",DK 
print "B   (kG) = ",Bkg
print "g        = ",g
kthz=0.0208366
JFM=JFMK*kthz
D=DK*kthz
mub=5.788*1e-2 #meV/T  1 T = 10 kG
mub=5.788*1e-3 #meV/ kG	
mub=5.788*1e-3/(0.086217) #KELVIN/ kG	
mub=5.788*1e-3*kthz/(0.086217) #THZ/ kG	
B=Bkg*mub*g
'''
B=0
JFM = 1.0
D = float(sys.argv[2])
num_magnons=[3]
for nmag in num_magnons:
	configs=[]
	if (nmag==0):
		config=N.ones((L),dtype=int)
		configs.append(config)
	if (nmag==1):
		for i in range(L):
			config=N.ones((L),dtype=int)
			config[i]=0
			configs.append(config)	
	if (nmag==2):
		for i in range(L):
			for j in range(i+1,L):
				config=N.ones((L),dtype=int)
				config[i]=0
				config[j]=0
				#print config
				configs.append(config)
		for i in range(L):
			config=N.ones((L),dtype=int)
			config[i]=-1
			#print config
			configs.append(config)	
	
	if (nmag==3):
		for i in range(L):
			for j in range(i+1,L):
				for k in range(j+1,L):
					config=N.ones((L),dtype=int)
					config[i]=0
					config[j]=0
					config[k]=0
					#print config
					configs.append(config)
		for i in range(L):
			for j in range(i+1,L):
				config=N.ones((L),dtype=int)
				config[i]=-1
				config[j]=0
				#print config
				configs.append(config)
		for i in range(L):
			for j in range(i+1,L):
				config=N.ones((L),dtype=int)
				config[i]=0
				config[j]=-1
				#print config
				configs.append(config)
	
	print "Number of configurations in # magnon sector",nmag," = ",len(configs)
	hilbert=len(configs)
	H=N.zeros((hilbert,hilbert),dtype=float)
	for i in range(hilbert):
		#print "setting up matrix elements for i = ",i
		diaghint=0.0
		for n in configs[i]:
			if (n==1): diaghint+=(-D - B)
			if (n==-1): diaghint+=(-D + B)
		for a in range(L):
			an=(a+1)%L
			diaghint+=(-JFM*(configs[i][a]*configs[i][an]))
			if ((configs[i][a]==1 and configs[i][an]==0) or (configs[i][a]==0 and configs[i][an]==1)): # 1 magnon exchange 
				newconfig=copy.deepcopy(configs[i])
				newconfig[a]=configs[i][an]
				newconfig[an]=configs[i][a]
				loc,present=search_config(newconfig,configs)
				if (present): H[i,loc]+=(-(JFM))

			if ((configs[i][a]==1 and configs[i][an]==-1) or (configs[i][a]==-1 and configs[i][an]==1)): # 1 -1 with 0 0 
				newconfig=copy.deepcopy(configs[i])
				newconfig[a]=0
				newconfig[an]=0
				loc,present=search_config(newconfig,configs)
				if (present): H[i,loc]+=(-(JFM))

			if ((configs[i][a]==0 and configs[i][an]==0)):     # 0 0 to 1 -1 and -1 1 
				newconfig=copy.deepcopy(configs[i])
				newconfig[a]=1
				newconfig[an]=-1
				loc,present=search_config(newconfig,configs)
				if (present): H[i,loc]+=(-(JFM))
				
				newconfig=copy.deepcopy(configs[i])
				newconfig[a]=-1
				newconfig[an]=1
				loc,present=search_config(newconfig,configs)
				if (present): H[i,loc]+=(-(JFM))

			#I started adding here
                        if ((configs[i][a]==0 and configs[i][an]==-1) or (configs[i][a]==-1 and configs[i][an]==0)): #0 -1 exchange
                                newconfig=copy.deepcopy(configs[i])
                                newconfig[a]=configs[i][an]
                                newconfig[an]=configs[i][a]
                                loc,present=search_config(newconfig,configs)
                                if (present): H[i,loc]+=(-(JFM))

		H[i,i]=diaghint

	eigs,vecs=N.linalg.eigh(H)
	print '3 magnon energy'
	for energy in eigs[:20]:
		print energy	
	correlator=three_magn_corr(configs,vecs[:,0])
        print "3 magnon correlator in the ground state"
        print "i", "j", "correlator"
	special=(L/2)-1
	aa=0
        for i in range(L):
                for j in range(L):
                       # if (i!=j) and (i!=special) and (j!=special):
			print i, j, correlator[i,j]
			aa+=1
	print "number of terms=", aa
	
	'''
	print " "
	
	print "3 magnon correlator in the 1st excited state"
	print "i", "j", "correlator"
	fcorr = three_magn_corr(configs,vecs[:,1])
	for i in range(L):
                for j in range(L):
                        print i, j, fcorr[i,j]
	
	print " "

	print "3 magnon correlator in the 3rd excited state"
        print "i", "j", "correlator"
        f3corr = three_magn_corr(configs,vecs[:,5])
        for i in range(L):
                for j in range(L):
                        print i, j, f3corr[i,j]
	'''
