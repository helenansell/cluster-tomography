'''
Created by: Helen Ansell
Last updated: 15 October 2023

Functions for cluster counting in cluster tomography

These functions are called by percolation_cluster_tomography.py
'''

import collections

# Get the cluster IDs at each site on a line through a percolation sample at angle tan(gamma) = 1/n
# x and y are the locations of one line endpoint
def get_angled_line_clusters(size, HK_Matrix, n, x, y=None):
	if n == 0:
		line_site_ids = [x+y*size  for y in range(size)]
	else:
		line_site_ids = [(a % size)+y*size for y in range(size) for a in range(x+y*n,x+(y+1)*n)]

	line_cluster_ids = [HK_Matrix[id] for id in line_site_ids]

	return line_cluster_ids
	
# Count number of clusters along a line segment
def count_clusters(line_cluster_ids):
	line_cluster_counts = collections.Counter(line_cluster_ids) # Number of sites in each cluster

	if 0 in line_cluster_counts: del line_cluster_counts[0] # Remove zero from the count
		
	cluster_count = len(line_cluster_counts) # Number of visited clusters
	
	return cluster_count


