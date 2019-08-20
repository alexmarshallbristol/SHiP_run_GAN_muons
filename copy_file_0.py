
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-jobid', action='store', dest='jobid', type=int,
					help='jobid')
results = parser.parse_args()
jobid = int(results.jobid)

list_of_file_ID = np.load('list_of_file_ID.npy')

print(jobid, list_of_file_ID[jobid])

file = '/eos/experiment/ship/user/amarshal/HUGE_GAN_random_id/huge_generation_april_%s.npy'%list_of_file_ID[jobid]

print(file)

command = 'cp %s muons.npy'%file

os.system(command)