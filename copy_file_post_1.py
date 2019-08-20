
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-jobid', action='store', dest='jobid', type=int,
					help='jobid')
results = parser.parse_args()
jobid = int(results.jobid)

jobid = jobid + 7500

list_of_file_ID = np.load('list_of_file_ID.npy')

print(jobid, list_of_file_ID[jobid])

file = '/eos/experiment/ship/user/amarshal/HUGE_GAN_random_id/huge_generation_april_%s.npy'%list_of_file_ID[jobid]

print(file)

command = 'cp ship.conical.MuonBack-TGeant4_rec.root /eos/experiment/ship/user/amarshal/HUGE_GAN_random_id_FairSHiP/GAN_ship.conical.MuonBack-TGeant4_rec_%s.root'%list_of_file_ID[jobid]

os.system(command)