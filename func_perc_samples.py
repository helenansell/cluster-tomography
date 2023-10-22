'''
Created by: Helen Ansell
Last updated: 15 October 2023

Functions for generating site and bond percolation samples on a 2d square lattice 
for a chosen linear size and occupation probability with free boundary conditions.

Functions for imposing periodic boundary conditions (pbc) in a chosen direction.

These functions are called by percolation_cluster_tomography.py

'''
import numpy as np

# Union and find are used as part of the HK routine to identify clusters
def find(x):
    global labels
    y = int(x)
    while (labels[y] != y):
        y = labels[y]
    while (labels[x] != x):
        z = labels[x]
        labels[x] = y
        x = z
    return y


def union(a, b):
    global labels
    labels[find(a)] = find(b)

# Generate a 2D site percolation sample
# Returns a matrix containing cluster IDs
# Size is linear size of sample, pbc_x and pbc_y are 1 if pbc are required in that direction and 0 otherwise
def get_2D_site_matrix(size, occ):
	global labels
	labels = [0]
	HK_Matrix = np.random.choice([0, 1], size=(size**2), p=[1 - occ, occ]).tolist()
	
	largest_label = 0

	for y in range(size):
		for x in range(size):
			site = y*size+x
		
			if x > 0:
				left = HK_Matrix[site-1]
			else :
				left = 0
				
			if y > 0:
				forward = HK_Matrix[site-size]
				
			else:
				forward = 0
			
			
			if HK_Matrix[site] == 1:

				if left != 0 and forward == 0:
					HK_Matrix[site] = find(left)
					

				elif left == 0 and forward != 0:
					HK_Matrix[site] = find(forward)
					
				elif left == 0 and forward == 0:
					largest_label += 1
					HK_Matrix[site] = largest_label
					labels.append(largest_label)
					
				elif left != 0 and forward != 0:
					union(left, forward)
					HK_Matrix[site] = find(left)
					

	for y in range(size):
		for x in range(size):
			site = y*size+x
			if HK_Matrix[site] != 0:
				HK_Matrix[site] = find(HK_Matrix[site])
											
	return HK_Matrix
    

# Generate a 2D bond percolation sample
# Returns a matrix containing cluster IDs and the left and forward bond matrices
# Size is linear size of sample, pbc_x and pbc_y are 1 if pbc are required in that direction and 0 otherwise
def get_2D_bond_matrix(size, occ):
	global labels
	labels = [0]
	largest_label = 0

	HK_Matrix = np.zeros(size*size).tolist()
	left_Matrix = np.random.choice([0, 1], size=(size**2), p=[1 - occ, occ]).tolist()
	forward_Matrix = np.random.choice([0, 1], size=(size**2), p=[1 - occ, occ]).tolist()
	
	
	for y in range(size):
		for x in range(size):
			#left, forward and up are positions in the left_matrix, forward_matrix and forward_matrix
			# hk_left, hk_forward and hk_up are positions in the HK_Matrix (which becomes the output)
			site = y*size+x
			
			if x > 0:
				left = left_Matrix[site-1]
				hk_left = HK_Matrix[site-1]
			
			else:
				left = 0
				hk_left = 0


			if y > 0:
				forward = forward_Matrix[site-size]
				hk_forward = HK_Matrix[site-size]
			
			else:
				forward = 0
				hk_forward = 0
				
			if left != 0 and forward == 0:
				HK_Matrix[site] = find(hk_left)

			elif left == 0 and forward != 0:
				HK_Matrix[site] = find(hk_forward)

			elif left == 0 and forward == 0:
				largest_label = largest_label + 1
				HK_Matrix[site] = largest_label
				labels = np.append(labels, largest_label)

			elif left != 0 and forward != 0:
				union(hk_left, hk_forward)
				HK_Matrix[site] = find(hk_left)
						
	
	for y in range(size):
		for x in range(size):
			site = y*size+x
			
			HK_Matrix[site] = find(HK_Matrix[site])

							
	return HK_Matrix, left_Matrix, forward_Matrix
	

# Apply pbc in the x direction for site percolation
def apply_pbc_x_site(size, HK_Matrix):
	
	for y in range(size):
		site = y*size
		left = HK_Matrix[site+(size-1)]
		HKsite = HK_Matrix[site]

		if left!= 0 and HKsite != 0:
			union(left, HKsite)

	for y in range(size):
		for x in range(size):
			site = y*size+x
			if HK_Matrix[site] != 0:
				HK_Matrix[site] = find(HK_Matrix[site])

# Apply pbc in the y direction for site percolation
def apply_pbc_y_site(size, HK_Matrix):

	for x in range(size):
		forward = HK_Matrix[(size-1)*size+x]			
		HKsite = HK_Matrix[x]
						
		if forward!= 0 and HKsite!= 0:
			union(forward, HKsite)

	for y in range(size):
		for x in range(size):
			site = y*size+x
			if HK_Matrix[site] != 0:
				HK_Matrix[site] = find(HK_Matrix[site])

# Apply pbc in the x direction for bond percolation
def apply_pbc_x_bond(size, HK_Matrix, left_Matrix):

	for y in range(size):
		site = y*size
		left = left_Matrix[site+size-1]
				
		if left != 0:
			union(HK_Matrix[site],HK_Matrix[site+size-1])

	for y in range(size):
		for x in range(size):
			site = y*size+x
			HK_Matrix[site] = find(HK_Matrix[site])

# Apply pbc in the y direction for bond percolation
def apply_pbc_y_bond(size, HK_Matrix, forward_Matrix):

	for x in range(size):
		site = (size-1)*size+x
		forward = forward_Matrix[site]
						
		if forward != 0:
			union(HK_Matrix[x],HK_Matrix[site])

	for y in range(size):
		for x in range(size):
			site = y*size+x
			HK_Matrix[site] = find(HK_Matrix[site])


# Apply pbc_x according to percolation type
def apply_pbc_x(type, size, HK_Matrix, left_Matrix=0):
	if type == "site":
		apply_pbc_x_site(size, HK_Matrix)

	elif type == "bond":
		apply_pbc_x_bond(size, HK_Matrix, left_Matrix)	

# Apply pbc_y according to percolation type
def apply_pbc_y(type, size, HK_Matrix, forward_Matrix=0):
	if type == "site":
		apply_pbc_y_site(size, HK_Matrix)

	elif type == "bond":
		apply_pbc_y_bond(size, HK_Matrix, forward_Matrix)			
