## THE GLOBAL CODE

import numpy as np
import os	
import argparse
import sys


def get_seq_pdb(filename,chain_name):

	V = list()
	count = 0
	
	A = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
	B = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

	with open(filename,"r") as f:
		for line in f:
			if line.startswith("ATOM") and "CA" in line:
				tokens = line.split()
				if chain_name == None:				
				
					if len(tokens[3]) == 3:
						V.append(tokens[3])

					if len(tokens[3]) == 4:
						if count == 0:				
							V.append(tokens[3][1:4])
							count = count +1
							continue
						if count == 1:
							count = 0						
				if chain_name != None:
					if chain_name == line[21]:				

						if len(tokens[3]) == 3:
							V.append(tokens[3])
					
						if len(tokens[3]) == 4:
							if count == 0:				
								V.append(tokens[3][1:4])
								count = count +1
								continue
							if count == 1:
								count = 0						
	
	
	if len(V) == 0:
		print("no residues found")
		quit()

	for i,x in enumerate(V):
		for j in range(20):		
			if (x in A):		
				if x == A[j]:
					V[i] = B[j] 						
			else:
				print("unknown amino acid residue(s) found")	
				quit()		
	
	W = ''.join(V)	
#	print(type(W))
#	print len(V), len(W)
	return W


def get_score(filename):

	"""
	Gets the score from the output of clustalw

	Parameters
	----------
	filename - score.txt in this case
	
	Returns
	-------
	S - Similarity Score 

	"""
	
	words = None
	with open(filename,"r") as f4:
		lines = f4.readlines()
		for line in lines:
			if "Sequences" in line and "Aligned." in line:
				words = line.split()
		score = float(words[4]) 
	return score


def get_coordinates_pdb(filename,chain_name):
	
	"""
	Gets coordinates as a two dimensional array from a given pdb file 
	Can get whole protein's coordinates as a single array or get any chain coordinates if specified		
	Considers duplicate atoms too!

	Parameters
	----------
	filename - pdb file

	Returns
	-------
	V : array
		(N,D) matrix, where N is points and D is dimension.
	
	"""
	x_column = None
	V = list()
	count =0

	with open(filename, 'r') as f:
		lines = f.readlines()
		for line in lines:
			if line.startswith("ATOM") and "CA" in line:
				tokens = line.split()
					
				if chain_name == None:
					if len(tokens[3]) == 3:
						try:
							x = line[30:38]
							y = line[38:46]
							z = line[46:54]
							V.append(np.asarray([x, y ,z], dtype=float))
						except:
							exit("error: Parsing input for the following line: \n{0:s}".format(line))
					
					if len(tokens[3]) == 4:
						count = count +1
	
						if count == 1:
							o1 = float(line[54:60])
							x1 = float(line[30:38])
							y1 = float(line[38:46])
							z1 = float(line[46:54])     
			
							continue
			
						if count ==2:
							o2 = float(line[54:60])
							x2 = float(line[30:38])
							y2 = float(line[38:46])
							z2 = float(line[46:54])       
			     				
							if o1 >= o2:
								X = x1  
								Y = y1 
								Z = z1
	
							if o1 < o2:
								X = x2 
								Y = y2
								Z = z2
							
							V.append(np.asarray((X,Y,Z), dtype = float))
							count = 0
						
				if chain_name != None:				
					if chain_name == line[21]:
						if len(tokens[3]) == 3:
							try:
								x = line[30:38]
								y = line[38:46]
								z = line[46:54]
								V.append(np.asarray([x, y ,z], dtype=float))
							except:
								exit("error: Parsing input for the following line: \n{0:s}".format(line))
						
						if len(tokens[3]) == 4:
							count = count +1
	
							if count == 1:
								o1 = float(line[54:60])
								x1 = float(line[30:38])
								y1 = float(line[38:46])
								z1 = float(line[46:54])     
								continue
				
							if count ==2:
								o2 = float(line[54:60])
								x2 = float(line[30:38])
								y2 = float(line[38:46])
								z2 = float(line[46:54])       
			     						
								if o1 >= o2:
									X = x1  
									Y = y1 
									Z = z1
		
								if o1 < o2:
									X = x2 
									Y = y2
									Z = z2
							
								V.append(np.asarray((X,Y,Z), dtype = float))
								count = 0									
						
	V = np.asarray(V)
	
	if len(V) == 0:
		print("No coordinates found")
		quit()
	
	#print len(V)	
	return V

def get_gap_array(filename,protein1,protein2):

	V = list()
	W = list()
	S = list()
	S1 = list()
	S2 = list()
	S = list()
	j = 1
	k = 1
	SA = list()
	SB = list()
	
	with open(filename, "r") as f1:
		for line in f1:	
			if line.startswith(protein1): 
				token1 = line.split()
				V.append(token1[1])
			if line.startswith(protein2):
				token2 = line.split()
				W.append(token2[1])
	
	sa = list(''.join(V))
	sb = list(''.join(W))
	#print len(sa),len(sb)
	

	for i,e in enumerate(sa):
		
		if e == "-":
			SA.append("-")
			continue
		else:
			SA.append(j)
			j = j+1

	for i,e in enumerate(sb):
		
		if e == "-":
			SB.append("-")
			continue
		else:
			SB.append(k)
			k = k+1

	if "-" in SA or "-" in SB:
	
		for i,e in enumerate(SA):
			if e == "-":	
				S1.append(i)
	
		for i in range(len(S1)):
			SB.pop(S1[i] - i)
			SA.pop(S1[i] - i)

	
		for i,e in enumerate(SB):
			if e == "-":	
				S2.append(i)
	
		for i in range(len(S2)):
			SA.pop(S2[i] - i)
			SB.pop(S2[i] - i)
	
	return SA,SB	
	
def centroid(X):
	
	"""
	Centroid is the mean position of all the points in all of the coordinate
	directions, from a vectorset X.

	C = sum(X)/len(X)

	Parameters
	----------
	X : array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	C : float
		centroid
	"""
	
	C = X.mean(axis=0)
	return C


def quaternion_rmsd(P, Q):
	
	"""
	Rotate matrix P unto Q and calculate the RMSD
	based on doi:10.1016/1049-9660(91)90036-O

	Parameters
	----------
	P : array
		(N,D) matrix, where N is points and D is dimension.
	Q : array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	rmsd : float
	"""

	rot = quaternion_rotate(P, Q)
	P = np.dot(P, rot)
	return rmsd(P, Q)


def quaternion_transform(r):

	"""
	Get optimal rotation
	note: translation will be zero when the centroids of each molecule are the
	same
	"""

	Wt_r = makeW(*r).T
	Q_r = makeQ(*r)
	rot = Wt_r.dot(Q_r)[:3, :3]
	return rot


def makeW(r1, r2, r3, r4=0):

	"""
	matrix involved in quaternion rotation
	"""

	W = np.asarray([
		[r4, r3, -r2, r1],
		[-r3, r4, r1, r2],
		[r2, -r1, r4, r3],
		[-r1, -r2, -r3, r4]])
	return W


def makeQ(r1, r2, r3, r4=0):

	"""
	matrix involved in quaternion rotation
	"""

	Q = np.asarray([
		[r4, -r3, r2, r1],
		[r3, r4, -r1, r2],
		[-r2, r1, r4, r3],
		[-r1, -r2, -r3, r4]])
	return Q


def quaternion_rotate(X, Y):

	"""
	Calculate the rotation

	Parameters
	----------
	X : array
		(N,D) matrix, where N is points and D is dimension.
	Y: array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	rot : matrix
		Rotation matrix (D,D)
	"""

	N = X.shape[0]
	W = np.asarray([makeW(*Y[k]) for k in range(N)])
	Q = np.asarray([makeQ(*X[k]) for k in range(N)])
	Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
	W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
	A = np.sum(Qt_dot_W, axis=0)
	eigen = np.linalg.eigh(A)
	r = eigen[1][:, eigen[0].argmax()]
	rot = quaternion_transform(r)
	return rot


def kabsch_rmsd(P, Q, translate=False):

	"""
	Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.

	Parameters
	----------
	P : array
		(N,D) matrix, where N is points and D is dimension.
	Q : array
		(N,D) matrix, where N is points and D is dimension.
	translate : bool
		Use centroids to translate vector P and Q unto each other.

	Returns
	-------
	rmsd : float
		root-mean squared deviation
	"""

	if translate:
		Q = Q - centroid(Q)
		P = P - centroid(P)

	P = kabsch_rotate(P, Q)
	return rmsd(P, Q)


def kabsch_rotate(P, Q):

	"""
	Rotate matrix P unto matrix Q using Kabsch algorithm.

	Parameters
	----------
	P : array
		(N,D) matrix, where N is points and D is dimension.
	Q : array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	P : array
		(N,D) matrix, where N is points and D is dimension,
		rotated

	"""
	U = kabsch(P, Q)

	# Rotate P
	P = np.dot(P, U)
	return P


def kabsch(P, Q):

	"""
	Using the Kabsch algorithm with two sets of paired point P and Q, centered
	around the centroid. Each vector set is represented as an NxD
	matrix, where D is the the dimension of the space.

	The algorithm works in three steps:
	- a centroid translation of P and Q (assumed done before this function
	  call)
	- the computation of a covariance matrix C
	- computation of the optimal rotation matrix U

	For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm

	Parameters
	----------
	P : array
		(N,D) matrix, where N is points and D is dimension.
	Q : array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	U : matrix
		Rotation matrix (D,D)
	"""

	# Computation of the covariance matrix
	C = np.dot(np.transpose(P), Q)

	# Computation of the optimal rotation matrix
	# This can be done using singular value decomposition (SVD)
	# Getting the sign of the det(V)*(W) to decide
	# whether we need to correct our rotation matrix to ensure a
	# right-handed coordinate system.
	# And finally calculating the optimal rotation matrix U
	# see http://en.wikipedia.org/wiki/Kabsch_algorithm
	V, S, W = np.linalg.svd(C)
	d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

	if d:
		S[-1] = -S[-1]
		V[:, -1] = -V[:, -1]

	# Create Rotation matrix U
	U = np.dot(V, W)

	return U


def rmsd(V, W):

	"""
	Calculate Root-mean-square deviation from two sets of vectors V and W.

	Parameters
	----------
	V : array
		(N,D) matrix, where N is points and D is dimension.
	W : array
		(N,D) matrix, where N is points and D is dimension.

	Returns
	-------
	rmsd : float
		Root-mean-square deviation between the two vectors
	"""

	D = len(V[0])
	N = len(V)
	result = 0.0
	for v, w in zip(V, W):
		result += sum([(v[i] - w[i])**2.0 for i in range(D)])
	return np.sqrt(result/N)


def takesecond(elem):
    return elem[1]


def main():

	if len(sys.argv) == 1:
		print "input format - python file.py 1HAB_A list.txt "
		quit()
	
	a = sys.argv[1]		
	x = sys.argv[2] 
	S = list()
	k = 0

	print "\n" "Proteins are sorted based on their RMSD values after superposition with", a

	with open(x,"r") as l:	
		for line in l:
			M = list()
			N = list()

			y = line.split()
			b = [a[:4].lower(),y[0][:4].lower()]		# b - list of two protein names in lower_case
		#	print b	
	
			if a[:4] != a[:4].upper() or y[0][:4] != y[0][:4].upper():
				print "protein name should be in caps - e.g. 1HAB_A "
				quit() 
		
			os.chdir("/home/guest/Downloads/sai/final_files/pdbs/")	  

			V = [b[0] + ".pdb", b[1] + ".pdb"] 								# for eg. V = ["1g4i.pdb", "1mku.pdb"]
	
			if len(a) == 4:	
				p_seq = get_seq_pdb(V[0],None)
				p_coords = get_coordinates_pdb(V[0],None) 
	
			elif len(a) == 6 and a[4] == "_":
				p_seq = get_seq_pdb(V[0],a[5])											
				p_coords = get_coordinates_pdb(V[0],a[5]) 

			else:
				print "\ninvalid format given in command line - ", sys.argv[1]
				print("Enter just protein name IN CAPS - (e.g) - 1HAB or 1HAB_A")
				quit()

			if len(y[0]) == 4:
				q_seq = get_seq_pdb(V[1],None)
				q_coords = get_coordinates_pdb(V[1],None)

			elif len(y[0]) == 6 and y[0][4] == "_":
				q_seq = get_seq_pdb(V[1],y[0][5])		
				q_coords = get_coordinates_pdb(V[1],y[0][5])

			else:
				print "\ninvalid format given in text file - ", sys.argv[2], y[0]
				print("Enter just protein name IN CAPS - (e.g) - 1HAB or 1HAB_A \n")
				quit()

			with open("super.txt","w") as f:
				f.write(">" + a)
				f.write("\n")	
				f.write(p_seq)
				f.write("\n")
				f.write(">" + y[0])
				f.write("\n")
				f.write(q_seq)
				f.write("\n") 	
	
			if a == y[0]:
				p_centroid = centroid(p_coords)
				q_centroid = centroid(q_coords)
		
				p_coords -= p_centroid
				q_coords -= q_centroid		
				
				rmsd_1 = round(quaternion_rmsd(p_coords,q_coords),4) 
				S.append((y[0], rmsd_1))

			else:
				os.system("clustalw2 super.txt > super_score.txt") 				#calling clustalw and storing its result in a text file 
					
				A,B = get_gap_array("super.aln",a,y[0])	
	
				for i in A:
					M.append(np.asarray(p_coords[i-1]))
 		
				for j in B:
					N.append(np.asarray(q_coords[j-1]))		

				M = np.asarray(M)
				N = np.asarray(N)
	
				p_centroid = centroid(M)
				q_centroid = centroid(N)
		
				M -= p_centroid
				N -= q_centroid
			
				if len(M) == len(N):
					rmsd_1 = round(quaternion_rmsd(M,N),4) 
					S.append((y[0], rmsd_1))
			
				else:
					print "fuck off, Couldn't superpose", y[0], "\n"
					continue		

		S.sort(key=takesecond)	
		
		
		print "\n"
		while k < len(S):
			print S[k][0], " RMSD - ", S[k][1]
			k = k + 1
		print "\n"
	

if __name__ == "__main__":
    main()



