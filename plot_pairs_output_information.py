
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
mpl.use('Agg')
import matplotlib.pyplot as plt
import math
from matplotlib.colors import LogNorm
plt.rc('text', usetex=True)
plt.rcParams['savefig.dpi'] = 100
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams.update({'font.size': 15})






'''
# pretend for now that 'collected_pair_info.npy' is the final merged output
'''







# [pair_weight, nmeas_i, nmeas_j, rchi2_i, rchi2_j, P_i, P_j, hits_before_and_after_i, hits_before_and_after_j, doca, fid, dist, xv, yv, zv, HNL_mom]

pair_information = np.load('collected_pair_info.npy')

print(np.shape(pair_information))

def get_errors(array, weights, start_value, end_value, bins):
	'''
	Return normalised plot with errors.
	'''
	bin_width = (end_value - start_value)/bins
	x = np.empty(bins)
	y = np.empty(bins)
	yerr = np.empty(bins)
	for i in range(0, bins):
		x[i] = i*bin_width + bin_width/2 + start_value
		to_delete = np.where(array > x[i] + bin_width/2)
		array_bin = np.delete(array,to_delete)
		weights_bin = np.delete(weights,to_delete)

		to_delete = np.where(array_bin < x[i] - bin_width/2)
		array_bin = np.delete(array_bin,to_delete)
		weights_bin = np.delete(weights_bin,to_delete)

		y[i] = np.sum(weights_bin)

		weights_bin_sq = weights_bin*weights_bin

		yerr[i] = math.sqrt(np.sum(weights_bin_sq))

	sum_y = np.sum(y)
	# print(sum_y)
	y = y/sum_y
	yerr = yerr/sum_y

	return x, y, yerr

def plot(x,y,w,name,xlabel):

	fig = plt.figure(figsize=(7,6))
	plt.errorbar(x, y, yerr=w,drawstyle='steps-mid',capsize=0, elinewidth=1.5, linewidth=1.5,label='Background')
	plt.xlabel(xlabel)
	plt.ylabel('a.u.')
	plt.ylim(0,1)
	# plt.yscale('log')
	plt.legend(loc='upper right', bbox_to_anchor=(1, 0.99), fontsize=15,framealpha=0,borderpad=0.1,labelspacing=0.1,handletextpad=0,borderaxespad=0,columnspacing=0,markerfirst=True)
	plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)
	plt.savefig('plots/%s'%name)
	plt.close('all')

def plot_log(track_momentum_gan,weights_gan,name,xlabel,min_value,max_value,num_bins):

	bins = np.logspace(np.log10(min_value),np.log10(max_value),num=num_bins)

	gan_hist_log = np.histogram(track_momentum_gan, weights = weights_gan, bins=bins, range=[0, int(max_value)])
	plot_gan_hist_log = np.empty((0,3))
	sum_weights = np.sum(weights_gan)
	for i in range(0, np.shape(gan_hist_log[0])[0]):
		sum_of_weights_squared = 0
		for x in range(0, np.shape(track_momentum_gan)[0]):
			if track_momentum_gan[x] > gan_hist_log[1][i] and track_momentum_gan[x] < gan_hist_log[1][i+1]:
				sum_of_weights_squared += weights_gan[x]**2
		error = math.sqrt(sum_of_weights_squared)
		plot_gan_hist_log = np.append(plot_gan_hist_log, [[(gan_hist_log[1][i]+gan_hist_log[1][i+1])/(2),gan_hist_log[0][i]/sum_weights,error/sum_weights]], axis=0)

	fig = plt.figure(figsize=(7,6))
	plt.errorbar(plot_gan_hist_log[:,0], plot_gan_hist_log[:,1], yerr=plot_gan_hist_log[:,2],drawstyle='steps-mid',capsize=0, elinewidth=1.5, linewidth=1.5,label='Background')
	plt.xlabel(xlabel)
	plt.ylabel('a.u.')
	plt.xlim(np.amin(track_momentum_gan),np.amax(track_momentum_gan))
	plt.ylim(0,1)
	plt.xscale('log')
	plt.legend(loc='upper right', bbox_to_anchor=(1, 0.99), fontsize=15,framealpha=0,borderpad=0.1,labelspacing=0.1,handletextpad=0,borderaxespad=0,columnspacing=0,markerfirst=True)
	plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)
	plt.savefig('plots/%s'%name)
	plt.close('all')

def plot_2d(x,y,w,name,xlabel,ylabel,bins,limits,log_z):

	cmap = plt.get_cmap('viridis')
	cmap.set_under(color='white')  

	fig = plt.figure(figsize=(7,6))
	if log_z == True:
		plt.hist2d(x, y,weights=w,bins=bins,range=limits,vmin=1E-7,cmap=cmap,norm=LogNorm())
	else:
		plt.hist2d(x, y,weights=w,bins=bins,range=limits,vmin=1E-7,cmap=cmap)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	cax = plt.axes([0.875, 0.15, 0.02, 0.7])
	cbar = plt.colorbar(cax=cax)
	plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)
	plt.savefig('plots/%s'%name)
	plt.close('all')



print('DOCA')
x, y, w = get_errors(pair_information[:,9], pair_information[:,0], 0, 800, 35)
plot(x, y, w,'Combi_DOCA','Distance of Closest Approach (cm)')

print('IP')
x, y, w = get_errors(pair_information[:,11], pair_information[:,0], 0, 1000, 35)
plot(x, y, w,'Combi_IP','Impact Parameter w.r.t Target (cm)')

print('vz')
x, y, w = get_errors(pair_information[:,14], pair_information[:,0], -10000, 10000, 35)
plot(x, y, w,'Combi_VZ','Z Coordinate of Reconstructed Vertex (cm)')




print('HNL mom')
x = pair_information[:,15]
w = pair_information[:,0]
plot_log(x, w, 'Combi_HNL_mom_log','Reconstructed Mother Momentum (GeV/c)',0.1,500.,35)




print('Vertex x:y')
x = pair_information[:,12]
y = pair_information[:,13]
w = pair_information[:,0]
plot_2d(x, y, w, 'Combi_Vertex_xy','Vertex X Coordinate (cm)','Vertex Y Coordinate (cm)',50,[[-1500,1500],[-1500,1500]],False)

print('Vertex z:y')
x = pair_information[:,14]
y = pair_information[:,13]
w = pair_information[:,0]
plot_2d(x, y, w, 'Combi_Vertex_zy','Vertex Z Coordinate (cm)','Vertex Y Coordinate (cm)',50,[[-5000,5000],[-1500,1500]],False)

















