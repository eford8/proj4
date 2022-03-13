#!/usr/bin/python3

from cmath import inf
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
		#score = random.random()*100;
		if(not banded): 
			n = len(seq1) + 1
			m = len(seq2) + 1
			score_array = np.empty([n, m], dtype = list)
			val = 0
			for i in range(0,n) :
				score_array[i,0] = [val, 0,0]
				val += 5
			val = 0
			for i in range(0,m) :
				score_array[0,i] = [i *5, 0,0]
				val += 5
			for i in range(1,n):
				for j in range(1,m):
					#print(score_array[i-1, j-1][0])
					sub = score_array[i-1, j-1][0] + SUB
					if seq1[i - 1] == seq2[j - 1] :
						sub = score_array[i-1, j-1][0] + MATCH
					indel_up = score_array[i-1, j][0] + INDEL
					indel_left = score_array[i, j - 1][0] + INDEL
					if(indel_left <= indel_up and indel_left <= sub):
						score_array[i,j] = [indel_left, i, j - 1]
					elif(indel_up <= indel_left and indel_up <= sub):
						score_array[i,j] = [indel_up, i - 1, j]
					else:
						score_array[i,j] = [sub, i-1, j-1]
			
			#print(score_array)
			score = score_array[n - 1,m - 1][0]
		else :
			n = len(seq2) + 1
			if (len(seq1) < len(seq2)) :
				n = len(seq1) + 1
			m = 7
			score_array = np.empty([n, m], dtype = list)
			for i in range(0,4) :
				#score_array[i,0] = [i * 5, 0,0]
				score_array[i,0] = i * 5
			for i in range(0,4) :
				#score_array[0,i] = [i * 5, 0,0]
				score_array[0,i] = i * 5
			for i in range(1,n):
				for j in range(0,m):
					#print(score_array[i-1, j-1][0])
					sub = inf
					indel_left = inf
					indel_up = inf
					if(i - 1 >= 0 and score_array[i-1, j]):
						#sub = score_array[i-1, j][0] + SUB
						sub = score_array[i-1, j] + SUB
						if seq1[i - 1] == seq2[j - 1] :
							#sub = score_array[i-1, j][0] + MATCH
							sub = score_array[i-1, j] + MATCH
					if(i - 1 >= 0 and j + 1 < m and score_array[i-1, j + 1]) :
						#indel_up = score_array[i-1, j+1][0] + INDEL
						indel_up = score_array[i-1, j+1] + INDEL
					if(j - 1 >= 0) :
						#indel_left = score_array[i, j - 1][0] + INDEL
						indel_left = score_array[i, j - 1] + INDEL

					if(i < 4 and j == 0) :
						continue
					elif(indel_left <= indel_up and indel_left <= sub):
						#score_array[i,j] = [indel_left, i, j - 1]
						score_array[i,j] = indel_left
					elif(indel_up <= indel_left and indel_up <= sub):
						#score_array[i,j] = [indel_up, i - 1, j+1]
						score_array[i,j] = indel_up
					else:
						#score_array[i,j] = [sub, i-1, j]
						score_array[i,j] = sub
			
			#print(score_array)
			#score = score_array[n - 1,m - 1][0]
			score = score_array[n - 1,m - 1]




		alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
			len(seq1), align_length, ',BANDED' if banded else '')
		alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
			len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
