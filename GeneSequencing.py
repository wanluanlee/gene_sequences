#!/usr/bin/python3

#from PyQt5.QtCore import QLineF, QPointF

import math

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

    def __init__( self ):
        pass

# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment
    def align( self, sequences, table, banded, align_length ):
        self.banded = banded
        self.MaxCharactersToAlign = align_length
        results = []

        for i in range(len(sequences)):
            jresults = []
            for j in range(len(sequences)):
                if j < i:
                   s = {}
                else:
                    if(banded == True):
                        if((abs(len(sequences[i][:align_length]) -
                                len(sequences[j][:align_length]))) > MAXINDELS):
                            score = float(math.inf)
                            alignment1 = 'No Alignment Found'
                            alignment2 = 'No Alignment Found'
                        else:
                            #myResult = self.bandit(sequences[7][:align_length], sequences[2][:align_length])
                            if(len(sequences[i][:align_length]) >= len(sequences[j][:align_length])):
                                myResult = self.bandit(sequences[j][:align_length], sequences[i][:align_length])
                            else:
                                myResult = self.bandit(sequences[i][:align_length], sequences[j][:align_length])

                            score = myResult[0]
                            alignment1 = myResult[1]
                            alignment2 = myResult[2]
                    else:
                        myResult = self.tofit(sequences[i][:align_length], sequences[j][:align_length])
                        score = myResult[0]
                        alignment1 = myResult[1]
                        alignment2 = myResult[2]

                    #self.tofit(sequences[0],sequences[1])
                    #self.bandit(sequences[0],sequences[1])
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
                    #score = i + j
                    #alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
                        #len(sequences[i]), align_length, ',BANDED' if banded else '')
                    #alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
                        #len(sequences[j]), align_length, ',BANDED' if banded else '')
###################################################################################################
                    s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
                    table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
                    table.repaint()
                jresults.append(s)
            results.append(jresults)
        return results

    def tofit(self, sequence1,sequence2):
        myMatrix = [[ 0 for i in range(len(sequence1) + 1) ] for j in range(len(sequence2) + 1)]
        # trace the path
        path = [[ 0 for i in range(len(sequence1) + 1) ] for j in range(len(sequence2) + 1)]
        for i in range(len(sequence1)+1):
            for j in range(len(sequence2)+1):
                # filled the top and the most left column
                myMatrix[0][i] = 5*i
                myMatrix[j][0] = 5*j
                # filled the top and the most left path as NULL
                path[0][i] = 'L'
                path[j][0] = 'T'
        path[0][0] = "NULL"
        for i in range(1,len(sequence2)+1):
            for j in range(1,len(sequence1)+1):
                # updating numbers
                corner = myMatrix[i-1][j-1]
                left = myMatrix[i][j-1]
                top = myMatrix[i-1][j]
                if(sequence1[j - 1] == sequence2[i - 1]):
                    comp1 = corner + MATCH
                else:
                    comp1 = corner + SUB
                comp2 = left + INDEL
                comp3 = top + INDEL
                minNum = min([comp1,comp2,comp3])
                myMatrix[i][j] = minNum
                if(comp1 == minNum):
                    path[i][j] = 'C'

                elif(comp2 == minNum):
                    path[i][j] = 'L'

                elif (comp3 == minNum):
                    path[i][j] = 'T'

        value = myMatrix[-1][-1]

        #Track back the path
        string_i = []
        string_j = []
        current = path[-1][-1]
        currentj = len(sequence1)
        currenti = len(sequence2)
        index_j = len(sequence1) -1
        index_i = len(sequence2) -1
        while(current != 0 and current != "NULL"):
            if (current == 'C'):
                string_j.append(sequence1[index_j])
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                index_j = index_j - 1
                currenti = currenti - 1
                currentj = currentj - 1
            if (current == 'T'):
                string_j.append('-')
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                currenti = currenti - 1
            if (current == 'L'):
                string_i.append(sequence1[index_i])
                string_j.append('-')
                index_j = index_j - 1
                currentj = currentj - 1

            current = path[currenti][currentj]
        string_i = string_i[::-1]
        string_j = string_j[::-1]

        string_i = string_i[:100]
        string_j = string_j[:100]

        return[value,string_i,string_j]

    def bandit(self,sequence1,sequence2):
        myMatrix = [[100 for i in range(7)] for j in range(len(sequence2) + 1)]
        path = [[100 for i in range(7)] for j in range(len(sequence2) + 1)]
        lenDiff = abs(len(sequence1) - len(sequence2))
        #initil first three
        for i in range (4):
            if(i == 0):
                myMatrix[i][0] = 'NULL'
                myMatrix[i][1] = 'NULL'
                myMatrix[i][2] = 'NULL'
                myMatrix[i][3] = 0
                myMatrix[i][4] = INDEL*1
                myMatrix[i][5] = INDEL*2
                myMatrix[i][6] = INDEL*3
                path[i][0] = 'NULL'
                path[i][1] = 'NULL'
                path[i][2] = 'NULL'
                path[i][3] = 0
                path[i][4] = 'L'
                path[i][5] = 'L'
                path[i][6] = 'L'

            if(i == 1):
                myMatrix[i][0] = 'NULL'
                myMatrix[i][1] = 'NULL'
                myMatrix[i][2] = INDEL*1

                path[i][0] = 'NULL'
                path[i][1] = 'NULL'
                path[i][2] = 'T'
            if(i == 2):
                myMatrix[i][0] = 'NULL'
                myMatrix[i][1] = INDEL*2

                path[i][0] = 'NULL'
                path[i][1] = 'T'
            if(i == 3):
                myMatrix[i][0] = INDEL*3
                path[i][0] = 'T'
            #----initaial buttom-----------
            for i in range(lenDiff + 1):
                myMatrix[len(sequence2) - i][6] = 'NULL'
                myMatrix[len(sequence2) - i][5] = 'NULL'
                myMatrix[len(sequence2) - i][4] = 'NULL'
            myMatrix[len(sequence2) - lenDiff - 1][6] = 'NULL'
            myMatrix[len(sequence2) - lenDiff - 1][5] = 'NULL'
            myMatrix[len(sequence2) - lenDiff - 2][6] = 'NULL'

        #Doing Operation
        index_i = 0
        index_j = 0
        toCompar = []
        toRight = 0
        ifNormal = False
        comp1 = "NULL"
        comp2 = "NULL"
        comp3 = "NULL"
        for i in range(len(sequence2) + 1):
            if(i > 3 and i < len(sequence2)):
                index_i = index_i + 1
                index_j = index_j + 1
            if(i == (len(sequence2)- lenDiff) + 1):
                ifNormal = True
            toRight = -1
            for j in range(7):
                current = myMatrix[i][j]
                if(current != "NULL"):
                    toRight = toRight + 1
                if(current == 100):
                    #check if edage case
                    if(ifNormal == False):
                        if((i-1) < 0):
                            top = "NULL"
                        else:
                            top = myMatrix[i - 1][j]
                            if (sequence1[index_j + toRight - 1] == sequence2[i - 1]):
                                comp1 = top + MATCH
                            else:
                                comp1 = top + SUB
                            toCompar.append(comp1)

                        if((j-1)<0):
                            left = "NULL"
                        else:
                            left = myMatrix[i][j - 1]
                            comp2 = left + INDEL
                            toCompar.append(comp2)
                        if((i-1) < 0 or (j+1) > 6):
                         corner = "NULL"
                        else:
                            corner = myMatrix[i - 1][j + 1]
                            comp3 = corner + INDEL
                            toCompar.append(comp3)
                    else:
                        if ((i - 1) < 0):
                            top = "NULL"
                        else:
                            top = myMatrix[i - 1][j]
                            comp1 = top + INDEL
                            toCompar.append(comp1)

                        if ((j - 1) < 0):
                            left = "NULL"
                        else:
                            left = myMatrix[i][j - 1]
                            comp2 = left + INDEL
                            toCompar.append(comp2)
                        if ((i - 1) < 0 or (j - 1) < 0):
                            corner = "NULL"
                        else:
                            corner = myMatrix[i - 1][j - 1]
                            if (sequence1[index_j + toRight - 1] == sequence2[i - 1]):
                                comp3 = corner + MATCH
                            else:
                                comp3 = corner + SUB
                            toCompar.append(comp3)

                    #-----------------------------------
                    minNum = min(toCompar)
                    myMatrix[i][j] = minNum
                    if (minNum == comp1 and ifNormal == False):
                        path[i][j] = "T"
                    elif (minNum == comp1 and ifNormal == True):
                        path[i][j] = "TN"
                    elif (minNum == comp3 and ifNormal == False):
                        path[i][j] = "C"
                    elif (minNum == comp3 and ifNormal == True):
                        path[i][j] = "CN"
                    elif (minNum == comp2 and ifNormal == False):
                        path[i][j] = "L"
                    elif (minNum == comp2 and ifNormal == True):
                        path[i][j] = "LN"

                    toCompar = []

        value = myMatrix[len(sequence2)][3]
        string_i = []
        string_j = []

        # Track back the path
        string_i = []
        string_j = []
        current = path[len(sequence2)][3]
        currentj = 3
        currenti = len(sequence2)
        index_j = len(sequence1) - 1
        index_i = len(sequence2) - 1
        while(current != 0 and current != "NULL"):
            if (current == 'T'):
                string_j.append(sequence1[index_j])
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                index_j = index_j - 1
                currenti = currenti - 1
            if (current == 'C'):
                string_j.append('-')
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                currentj = currentj + 1
                currenti = currenti - 1
            if (current == 'L'):
                string_j.append(sequence1[index_j])
                string_i.append('-')
                index_j-=1
                currentj = currentj - 1
            if(current == 'CN'):
                string_j.append(sequence1[index_j])
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                index_j = index_j - 1
                currenti = currenti - 1
                currentj = currentj - 1
            if(current == 'TN'):
                string_j.append('-')
                string_i.append(sequence2[index_i])
                index_i = index_i - 1
                currenti = currenti - 1
            if(current == 'LN'):
                string_j.append(sequence1[index_j])
                string_i.append('-')
                index_j = index_j - 1
                currentj = currentj - 1
            #print(path[currenti][currentj],"i= ",currenti, "j= ",currentj)
            current = path[currenti][currentj]


        string_i = string_i[::-1]
        string_j = string_j[::-1]

        string_i = string_i[:100]
        string_j = string_j[:100]

        return [value, string_i, string_j]



