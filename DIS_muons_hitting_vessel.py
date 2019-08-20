from __future__ import print_function
# example for accessing smeared hits and fitted tracks
import ROOT,os,sys,getopt
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
from rootpyPickler import Pickler
from decorators import *
import shipRoot_conf
shipRoot_conf.configure()
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import math
import glob 
import pickle 
import argparse
import shutil

files = glob.glob('/eos/experiment/ship/user/amarshal/HUGE_GAN_random_id_FairSHiP/GAN_ship.conical.MuonBack-TGeant4_rec_*')

print(np.shape(files))

hits = np.empty((0,6))
file_counter = 0

np.save('DIS_vessel_hits',hits)


hist_total = np.histogram2d(np.sqrt(np.add(hits[:,3]**2,np.add(hits[:,4]**2,hits[:,5]**2))),np.sqrt(np.add(hits[:,3]**2,hits[:,4]**2)),bins=50,range=[[0,350],[0,12]])
total = np.zeros((50,50))
# print(np.shape(hist_total[0]))
# quit()
for file in files:
	hits = np.empty((0,6))
	file_counter += 1
	if file_counter % 100 == 0: print(file_counter)
	f = ROOT.TFile(file,"read")
	tree = f.cbmsim
	N = tree.GetEntries()
	# print(file, N)
	for i in xrange(N):
		# if i % 1000 == 0: print(i, '/',N)
		tree.GetEntry(i)

		vessel_wall_hits = tree.vetoPoint

		for e in vessel_wall_hits:
			hits = np.append(hits,[[e.GetX(), e.GetY(), e.GetZ(), e.GetPx(), e.GetPy(), e.GetPz()]],axis=0)
			break

	hist = np.histogram2d(np.sqrt(np.add(hits[:,3]**2,np.add(hits[:,4]**2,hits[:,5]**2))),np.sqrt(np.add(hits[:,3]**2,hits[:,4]**2)),bins=50,range=[[0,350],[0,12]])
	# print(np.shape(hist[0]))
	total = np.add(hist[0],total)


	if file_counter % 500 == 0:

		plt.figure(figsize=(6,6))
		plt.imshow(np.flipud(total.T),norm=LogNorm(),extent=[0,350,0,12],aspect='auto',interpolation='none')
		plt.colorbar()
		plt.savefig('DIS')
		print(np.sum(total))
	# 	current = np.load('DIS_vessel_hits.npy')
	# 	hits = np.concatenate((current,hits),axis=0)
	# 	np.save('DIS_vessel_hits',hits)
	# 	hits = np.empty((0,6))

	# if file_counter % 1000 == 0: 
	# 	plt.hist2d(np.sqrt(np.add(hits[:,3]**2,np.add(hits[:,4]**2,hits[:,5]**2))),np.sqrt(np.add(hits[:,3]**2,hits[:,4]**2)),bins=50,norm=LogNorm(),range=[[0,350],[0,12]])
	# 	plt.colorbar()
	# 	plt.savefig('DIS')


# print(np.shape(hits))
np.save('DIS_vessel_hits',total)

# plt.hist2d(np.sqrt(np.add(hits[:,3]**2,np.add(hits[:,4]**2,hits[:,5]**2))),np.sqrt(np.add(hits[:,3]**2,hits[:,4]**2)),bins=50,norm=LogNorm(),range=[[0,350],[0,12]])
# plt.colorbar()
# plt.savefig('DIS')



print(np.sum(total))
plt.figure(figsize=(6,6))
plt.imshow(np.flipud(total.T),norm=LogNorm(),extent=[0,350,0,12],aspect='auto',interpolation='none')
plt.colorbar()
plt.savefig('DIS')
