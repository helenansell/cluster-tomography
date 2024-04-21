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

# Generate a 2D site percolation sample on a square lattice with free boundary conditions.
# Size is linear size of sample, occ is the occupancy probability for each site. 
#
# Returns a 1D matrix where each position in the matrix corresponds to a lattice site, indexed sequentially
# from top left to bottom right. Positions 0:size-1 correspond to the top row of the sample, positions
# 0,size,2*size...size*(size-1) correspond to the first column, etc. A value of 0 at a given location means the site
# is unoccupied, while values > 0 indicate the ID of the cluster the lattice site belongs to.
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
    

# Generate a 2D bond percolation sample on a square lattice with free boundary conditions.
# Returns HK_Matrix, left_Matrix, forward_Matrix.
# Size is the linear size of the sample, occ is the occupancy probability for each bond.
#
# HK_Matrix gives information about the clusters in the sample. Each position in the matrix corresponds to a 
# lattice site, indexed sequentially from top left to bottom right. Positions 0:size-1 correspond to the top row 
# of the sample, positions 0,size,2*size...size*(size-1) correspond to the first column, etc. The value at each 
# position gives the ID of the cluster the lattice site belongs to. 
# left_Matrix and forward_Matrix contain information about the connectivity of each lattice site. The sites are 
# indexed in the same way as HK_Matrix. A value of 1 in left_Matrix indicates that the bond between a site (with id 
# i) and its left neighbor (with id i-1) is present, while 0 means the bond is unoccupied. Similarly, in 
# forward_Matrix a value 1 means there is a bond between sites i and i-size, while 0 indicates there is no bond.
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
