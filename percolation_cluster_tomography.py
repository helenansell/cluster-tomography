'''
Created by: Helen Ansell
Last updated: 21 April 2024

Requires: func_perc_samples.py and func_line_clusters.py

Calculate the difference in cluster number count for the chosen tomography line type 
(as defined in fig. 1 of Ansell, Frank and Kovacs, Phys. Rev. Res. 5, 043218 (2023))
for site or bond percolation on a 2d square lattice.

Returns the mean value of DeltaN, the contribution to the cluster number count after the area law has been canceled.
For line types with angle dependence, results are calculated for angles gamma = arctan(1/n) for n <= 20

Note that to extract the value of b for line types 7, 8 and 10, contributions from lines of types 6 and 9 
must be subtracted from the values returned here. These additional lines are required so that the contour 
forms a closed loop when the boundary changing procedure is applied (see fig. S1 of the manuscript).
'''

import numpy as np
from time import time
import func_perc_samples as perc
import func_line_clusters as lc

# Parameters
perc_type = "site" #percolation type - must be "site" or "bond"
size = 32 #Linear side length of a sample
line_type=7 #Line types 1-10 
n_samples = 5 #Number of independent samples to analyze



half_size = int(size/2)

deltaNvals = []

start = time()

print("Performing cluster tomography on %i %s percolation samples of size %i for line type %i" %(size, perc_type, n_samples, line_type))

for sampleNo in range(n_samples):

	# Create a sample with free boundary conditions all around
	if perc_type == "site": 
		occ = 0.59275 # This is the occupancy probability, here set to the critical value for 2d site percolation
		HK_Matrix, left_Matrix, forward_Matrix = perc.get_2D_site_matrix(size, occ), 0, 0
		

	elif perc_type == "bond": 
		occ = 0.5000 # This is the occupancy probability, here set to the critical value for 2d bond percolation
		HK_Matrix, left_Matrix, forward_Matrix = perc.get_2D_bond_matrix(size, occ)

	else: 
		print("Perc_type must be 'site' or 'bond'")
		break

	# Calculate Delta N for the chosen line type
	

	# Line 1 - both ends in bulk
	# Note that for better statistics in a given sample, the two half-lines can be placed at
	# any positions in the sample where together the two form a closed contour 
	# Here, for simplicity, we demonstrate a subset of such posiitons
	if line_type == 1:
		# Apply pbc in x and y directions
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)

		line_1_counts = []

		for y in range(size):
			line_1a_clusters = HK_Matrix[y*size:y*size+half_size]
			line_1b_clusters = HK_Matrix[y*size+half_size:(y+1)*size]
			line_1A_clusters = HK_Matrix[y*size:(y+1)*size]

			line_1a_count = lc.count_clusters(line_1a_clusters)
			line_1b_count = lc.count_clusters(line_1b_clusters)
			line_1A_count = lc.count_clusters(line_1A_clusters)

			line_1_counts.append(line_1a_count + line_1b_count - line_1A_count)

		deltaNvals.append(np.mean(line_1_counts)/2)

	# Line 2 - partial surface line both ends on edge
	# Note that for better statistics in a given sample, the endpoints of the two half-lines 
	# can be placed at each of L positions top and bottom
	# Here we demonstrate a single placement top and bottom for simplicity
	elif line_type == 2:
		# Apply pbc in x direction 
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)

		# Count number of clusters along half-lines along top and bottom
		line_2a_clusters = HK_Matrix[:half_size]
		line_2b_clusters = HK_Matrix[half_size:size]
		line_2c_clusters = HK_Matrix[-size:-half_size]
		line_2d_clusters = HK_Matrix[-half_size:]

		line_2a_count = lc.count_clusters(line_2a_clusters)
		line_2b_count = lc.count_clusters(line_2b_clusters)
		line_2c_count = lc.count_clusters(line_2c_clusters)
		line_2d_count = lc.count_clusters(line_2d_clusters)

		# Count clusters along full lines along top and bottom
		line_2A_clusters = HK_Matrix[:size]
		line_2B_clusters = HK_Matrix[-size:]

		line_2A_count = lc.count_clusters(line_2A_clusters)
		line_2B_count = lc.count_clusters(line_2B_clusters)

		# delta N is half the difference between the sum of the two half counts and the full counts
		deltaNvals.append(np.mean([line_2a_count + line_2b_count - line_2A_count , line_2c_count + line_2d_count - line_2B_count])/2)

	# Line 3 - full surface line
	elif line_type ==  3:
		# Count number of clusters along top and bottom
		line_3a_clusters = HK_Matrix[:size]
		line_3b_clusters = HK_Matrix[-size:]

		line_3a_count = lc.count_clusters(line_3a_clusters)
		line_3b_count = lc.count_clusters(line_3b_clusters)

		# Apply pbc in x direction
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)

		# Count clusters again
		line_3a_clusters = HK_Matrix[:size]
		line_3b_clusters = HK_Matrix[-size:]

		line_3a_count -= lc.count_clusters(line_3a_clusters)
		line_3b_count -= lc.count_clusters(line_3b_clusters)

		deltaNvals.append(np.mean([line_3a_count, line_3b_count]))

	# Line 4 - one end in corner
	elif line_type ==  4:
		# Count number of clusters along half-lines along top and bottom
		line_4a_clusters = HK_Matrix[:half_size]
		line_4b_clusters = HK_Matrix[half_size:size]
		line_4c_clusters = HK_Matrix[-size:-half_size]
		line_4d_clusters = HK_Matrix[-half_size:]

		line_4a_count = lc.count_clusters(line_4a_clusters)
		line_4b_count = lc.count_clusters(line_4b_clusters)
		line_4c_count = lc.count_clusters(line_4c_clusters)
		line_4d_count = lc.count_clusters(line_4d_clusters)

		# Apply pbc in x direction
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)

		# Count clusters again
		line_4A_clusters = HK_Matrix[:size]
		line_4B_clusters = HK_Matrix[-size:]

		line_4A_count = lc.count_clusters(line_4A_clusters)
		line_4B_count = lc.count_clusters(line_4B_clusters)

		# delta N is half the difference between the sum of the two half counts and the full counts after applying pbc
		deltaNvals.append(np.mean([line_4a_count + line_4b_count - line_4A_count , line_4c_count + line_4d_count - line_4B_count])/2)

	# Line 5 - diagonal line
	elif line_type ==  5:
		# Count number of clusters along the two diagonal lines in each sample
		line_5A_clusters = [HK_Matrix[site] for site in range(0,size*size,size+1)]
		line_5B_clusters = [HK_Matrix[site] for site in range(size-1,size*(size-1)+1,size-1)]

		line_5A_count = lc.count_clusters(line_5A_clusters)
		line_5B_count = lc.count_clusters(line_5B_clusters)

		# Apply pbc in both directions
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)

		# Count number of clusters after bc change
		line_5A_clusters = [HK_Matrix[site] for site in range(0,size*size,size+1)]
		line_5B_clusters = [HK_Matrix[site] for site in range(size-1,size*(size-1),size-1)]

		line_5A_count -= lc.count_clusters(line_5A_clusters)
		line_5B_count -= lc.count_clusters(line_5B_clusters)

		deltaNvals.append(np.mean([line_5A_count,line_5B_count]))

	# Line 6 - line through bulk touching opposite edges
	elif line_type == 6:
		# Apply pbc in x direction
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)

		line_6_counts = [[] for _ in range(21)]
		for n in range(21):
			# start the line at each possible location
			line_6_x_counts = []

			for x in range(size):
				line_6_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, x)
				line_6_x_counts.append(lc.count_clusters(line_6_clusters))
			
			line_6_counts[n] = line_6_x_counts

		# Apply pbc in y direction 
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)

		for n in range(21):
			for x in range(size):
				line_6_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, x)
				line_6_counts[n][x] -= lc.count_clusters(line_6_clusters)

		deltaNvals.append([np.mean(a) for a in line_6_counts])


	# Line 7 - line through bulk touch adjacent edges
	elif line_type == 7:
		line_7_counts = [0 for _ in range(1,21)]

		for n in range(1,21):
			line_7_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, half_size)

			line_7_counts[n-1] = lc.count_clusters(line_7_clusters[:half_size]) + lc.count_clusters(line_7_clusters[-half_size:])

			if n > 1:
				for m in range(n-1):
					line_7_stripe_clusters = line_7_clusters[m*size+half_size:(m+1)*size+half_size]
					line_7_counts[n-1] += lc.count_clusters(line_7_stripe_clusters)

		# Apply pbc in x and y directions
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)

		for n in range(1,21):
			line_7_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, half_size)
			line_7_counts[n-1] -= lc.count_clusters(line_7_clusters)
			

		deltaNvals.append([line_7_counts[n]/2 for n in range(20)])
	
	# Line 8 - line through bulk from corner to edge
	# Note that to get the line 8 corner contribution, the contribution from internal stripes (line type 6) must be subtracted
	# The code here does not do this cancellation
	elif line_type == 8:
		line_8_counts = [0 for _ in range(2,21)]

		for n in range(2,21):
		# For each angle, only the first and last 'stripe' of line are of interest to us
		# The rest are just variants of standard full tomography
			line_8_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, 0)

			for m in range(n):
				line_8_stripe_clusters = line_8_clusters[m*size:(m+1)*size]
				line_8_counts[n-2] += lc.count_clusters(line_8_stripe_clusters)

		# Apply pbc in x and y directions
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)	

		for n in range(2,21):
		# For each angle, only the first and last 'stripe' of line are of interest to us
		# The rest are just variants of standard full tomography
			line_8_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, 0)
			line_8_counts[n-2] -= lc.count_clusters(line_8_clusters)

		deltaNvals.append([line_8_counts[n]/2 for n in range(19)])


	# Line 9 - line from edge ending in the bulk
	elif line_type == 9:
		# Apply pbc in x direction
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)

		line_9_counts = [[] for _ in range(21)]

		for n in range(21):
			# start the line at each possible location
			line_9_x_counts = []

			for x in range(size):
				line_9_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, x)

				if n > 0:
					line_9_x_counts.append(lc.count_clusters(line_9_clusters[:n*half_size])+lc.count_clusters(line_9_clusters[n*half_size:]))
				else:
					line_9_x_counts.append(lc.count_clusters(line_9_clusters[:half_size])+lc.count_clusters(line_9_clusters[half_size:]))
							
			line_9_counts[n] = line_9_x_counts

		# Apply pbc in y direction 
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)

		for n in range(21):
			for x in range(size):
				line_9_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, x)
				line_9_counts[n][x] -= lc.count_clusters(line_9_clusters)

		deltaNvals.append([np.mean(line_9_counts[n])/2 for n in range(21)])
		

	# Line 10 - line from corner ending in the bulk
	# Note that to get the line 10 corner contribution, the contribution from internal stripes (line type 6 and 9) must be subtracted
	# The code here does not do this cancellation
	elif line_type == 10:
		line_10_counts = [0 for _ in range(1,21)]

		for n in range(1,21):
		# For each angle, only the first and last 'stripe' of line are of interest to us
		# The rest are just variants of standard full tomography
			line_10_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, 0)

			# Count along half-line from each corner and the other half-line
			line_10a_clusters = line_10_clusters[:half_size]
			line_10b_clusters = line_10_clusters[half_size:size]

			line_10_counts[n-1] = lc.count_clusters(line_10a_clusters) + lc.count_clusters(line_10b_clusters) 

			if n > 1:

				line_10c_clusters = line_10_clusters[-size:-half_size]
				line_10d_clusters = line_10_clusters[-half_size:]

				line_10_counts[n-1] += lc.count_clusters(line_10c_clusters) + lc.count_clusters(line_10d_clusters)

			if n > 2: # There are internal lines of type 6 that must be counted
				for m in range(n-2):
					line_10_stripe_clusters = line_10_clusters[(m+1)*size:(m+2)*size]
					line_10_counts[n-1] += lc.count_clusters(line_10_stripe_clusters)

		# Apply pbc in x and y directions
		perc.apply_pbc_x(perc_type, size, HK_Matrix, left_Matrix)
		perc.apply_pbc_y(perc_type, size, HK_Matrix, forward_Matrix)	

		for n in range(1,21):
		# For each angle, only the first and last 'stripe' of line are of interest to us
		# The rest are just variants of standard full tomography
			line_10_clusters = lc.get_angled_line_clusters(size, HK_Matrix, n, 0)
			line_10_counts[n-1] -= lc.count_clusters(line_10_clusters)

		deltaNvals.append([line_10_counts[n]/2 for n in range(20)])


print("\nCalculation complete\n")

# Calculate the mean <Delta N> and standard error on the mean
if line_type < 6:
	meanDeltaN = np.mean(deltaNvals)
	steDeltaN = np.std(deltaNvals)/np.sqrt(n_samples)

	print("<Delta N > = %f ± %f" %(meanDeltaN,steDeltaN))
	

else:
	nAngles = len(deltaNvals[0])

	meanDeltaN = [np.mean([a[i] for a in deltaNvals]) for i in range(nAngles)]
	steDeltaN = [np.std([a[i] for a in deltaNvals])/np.sqrt(n_samples) for i in range(nAngles)]

	print("n \t <Delta N >")

	for n in range(21-nAngles, 21):
		print("%i \t %f ± %f" %(n, meanDeltaN[n-(21-nAngles)], steDeltaN[n-(21-nAngles)]))

end = time()

print("time elapsed = %f" %(end-start))

