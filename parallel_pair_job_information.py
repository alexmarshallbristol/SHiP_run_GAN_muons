import numpy as np
import math

track_location_array = np.load('track_location_array.npy')

number_of_fittracks = float(np.shape(track_location_array)[0])
# number_of_fittracks = 213.
number_of_jobs = 7500









print('Number of tracks:',int(number_of_fittracks))

print('number of pairs:',int((number_of_fittracks*(number_of_fittracks-1))/(2)))

print('with',number_of_jobs,'jobs, there will be',(number_of_fittracks*(number_of_fittracks-1))/(2*number_of_jobs),'pairs created in each.')

if (number_of_fittracks*(number_of_fittracks-1))/(2*number_of_jobs)-int((number_of_fittracks*(number_of_fittracks-1))/(2*number_of_jobs)) != 0:
	number_created_per = int(math.ceil((number_of_fittracks*(number_of_fittracks-1))/(2*number_of_jobs)))
else:
	number_created_per = int((number_of_fittracks*(number_of_fittracks-1))/(2*number_of_jobs))

print('So on all but the last job there will be',number_created_per,'pairs.')
print('and on the last there will be',int((number_of_fittracks*(number_of_fittracks-1))/(2)) - (number_created_per*(number_of_jobs-1)),'pairs.')

print('pair = [major, minor]')

total_pairs_created = 0


list_of_majors = range(0, int(number_of_fittracks))
number_of_minors = np.empty(np.shape(list_of_majors))
# print(np.shape(number_of_minors))
subtract_minors = 0
for i in range(0, np.shape(number_of_minors)[0]):
	number_of_minors[i] = int(np.shape(number_of_minors)[0] - 1 - subtract_minors) # -1 cause cannot have a pair with yourself
	subtract_minors += 1

overlap = 0
last_touched_major = 0
complete = False

job_order = np.zeros((int(number_of_jobs), 3))
#[job_id, start_major, start_minor]
job_order[0][0] = 0
job_order[0][1] = 0
job_order[0][2] = 1

for job_id in range(0,int(number_of_jobs)):

	if job_id != (number_of_jobs-1):
		total_pairs_created += number_created_per
	else:
		total_pairs_created = int((number_of_fittracks*(number_of_fittracks-1))/(2)) #all

	created_total = 0 

	for i in range(last_touched_major, np.shape(list_of_majors)[0]-1):
		if i == last_touched_major and overlap > 0:
			created_total += overlap
		else:
			created_total += number_of_minors[i]

		if total_pairs_created == int((number_of_fittracks*(number_of_fittracks-1))/(2)):
			complete = True
			break
		elif created_total > number_created_per:
			overlap = created_total - number_created_per
			
			last_touched_major = i
			break

	if complete == False:

		if overlap > 0:
			job_order[job_id+1][0] = job_id+1
			job_order[job_id+1][1] = last_touched_major
			job_order[job_id+1][2] = int(number_of_fittracks)-overlap

		else:
			job_order[job_id+1][0] = job_id+1
			job_order[job_id+1][1] = last_touched_major+1
			job_order[job_id+1][2] = 0

print(np.shape(job_order))
np.save('job_order',job_order)

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
# # OPEN up pretending to be on each job

# pairs_run_through = 0

# for job_id in range(0,int(number_of_jobs)):
# 	#pretend to be on job: job_id

# 	# print(' ')
# 	# print('On job',job_id)

# 	# print('job_order start:',job_order[job_id][1], job_order[job_id][2])
# 	try:
# 		# print('job_order end (not inclusive):',job_order[job_id+1][1], job_order[job_id+1][2])
# 		finish_major = job_order[job_id+1][1]
# 		finish_minor = job_order[job_id+1][2]
# 	except:
# 		# print('job_order end: GO TO END')
# 		finish_major = -1
# 		finish_minor = -1

# 	first_i = True

# 	for i in range(int(job_order[job_id][1]), int(number_of_fittracks)):

# 		if first_i == True:
# 			for j in range(int(job_order[job_id][2]), int(number_of_fittracks)):

# 				if i == finish_major and j == finish_minor:
# 					break

# 				else:
# 					pairs_run_through += 1
# 					# print([i,j])
# 		else:
# 			for j in range(i+1, int(number_of_fittracks)):

# 				if i == finish_major and j == finish_minor:
# 					break

# 				else:
# 					pairs_run_through += 1
# 					# print([i,j])

# 		first_i = False

# 		if i == finish_major and j == finish_minor:
# 				break



# print(pairs_run_through)













