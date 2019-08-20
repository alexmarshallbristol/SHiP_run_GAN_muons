
import numpy as np
import argparse
import os
# import os.path

parser = argparse.ArgumentParser()
parser.add_argument('-jobid', action='store', dest='jobid', type=int,
					help='jobid')
results = parser.parse_args()
jobid = int(results.jobid)

list_of_file_ID = np.load('list_of_file_ID.npy')

print(jobid, list_of_file_ID[jobid])


#check file exists

file_exists = os.path.exists('/eos/experiment/ship/user/amarshal/HUGE_GAN_random_id_FairSHiP/GAN_ship.conical.MuonBack-TGeant4_rec_%s.root'%list_of_file_ID[jobid])

if file_exists == True:
	file = '/eos/experiment/ship/user/amarshal/HUGE_GAN_random_id/huge_generation_april_%s.npy'%list_of_file_ID[jobid]

	print(file)

	command = 'cp %s muons.npy'%file

	os.system(command)
	
	f= open("file_is_not_present.txt","w+")
	print('File is not present... continue with job')
else:
	quit()