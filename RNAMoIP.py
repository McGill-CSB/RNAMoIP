###############################################################################
#Copyright (c) 2012, Vladimir Reinharz & Jerome Waldispuhl                    #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#* Redistributions of source code must retain the above copyright             #
#notice, this list of conditions and the following disclaimer.                #
#* Redistributions in binary form must reproduce the above copyright          #
#notice, this list of conditions and the following disclaimer in the          #
#documentation and/or other materials provided with the distribution.         #
#* Neither the name of the <organization> nor the                             #
#names of its contributors may be used to endorse or promote products         #
#derived from this software without specific prior written permission.        #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND#
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED#
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE       #
#DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY           #
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT   #
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS#
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE                  #
###############################################################################
from __future__ import with_statement
'''
@author: Vladimir Reinharz
   This program need gurobi installed.
   The input format is:
        1- The rna sequence
        2- The rna secondary structure (no pseudoknots)
        3- the path to the .desc files
        4- the fraction of basepairs that can be removed (in (0, 1))
        5- the max nb of components in motifs (Implemented for 4)

    It will output in sequence the first optimal solutions with the gurobi 
        general output.
    
'''

from gurobipy import *
import os
import sys
import re

NUMBER_THREADS = 0


def getArrayMotSequence(pos_let):
    """
        Should take as argument the second line of a .desc file
        This will create an array of strings, representing the different sequences of the motifs, in order
    """
    pos_let = pos_let[8:].split() #the first 10 characters are useless
    pos_let = [i.split('_') for i in pos_let] #now we have an array of tuples, in the first place the numb, and char in second
    if int(pos_let[0][0]) == 0:
        return None

    sequence = [str(pos_let[0][1])] #the firt sequence of the motif start with the first char

    i = 1
    while i < len(pos_let):
        if int(pos_let[i][0]) == 0:
            return None

        if int(pos_let[i][0]) - int(pos_let[i - 1][0]) > 5: #if the position does not follow, we start new sequence
            sequence.append(pos_let[i][1])
        elif int(pos_let[i][0]) - int(pos_let[i - 1][0]) == 1:
            sequence[len(sequence) - 1] = sequence[len(sequence) - 1] + pos_let[i][1] #if the position follow, we are on same sequence
        elif int(pos_let[i][0]) - int(pos_let[i - 1][0]) == 2:
            sequence[len(sequence) - 1] = sequence[len(sequence) - 1] + '.' + pos_let[i][1]
        elif int(pos_let[i][0]) - int(pos_let[i - 1][0]) == 3:
            sequence[len(sequence) - 1] = sequence[len(sequence) - 1] + '..' + pos_let[i][1]
        elif int(pos_let[i][0]) - int(pos_let[i - 1][0]) == 4:
            sequence[len(sequence) - 1] = sequence[len(sequence) - 1] + '...' + pos_let[i][1]
        elif int(pos_let[i][0]) - int(pos_let[i - 1][0]) == 5:
            sequence[len(sequence) - 1] = sequence[len(sequence) - 1] + '....' + pos_let[i][1]

        i += 1
    #Make sure at least 2/3 known
    not_known = 0
    tot = 0
    for i in sequence:
        if len(i) == 1:
            tot += 3
            not_known += 2
        elif len(i) == 2:
            tot += 3
            not_known += 1
        else:
            not_known = len(i.split('.')) - 1
            tot += len(i)
    if float(not_known) / tot >= 1.0 / 3 : 
        return None
    return sequence

def mergePositionsLists(List1, List2):
    """
        We merge the two arrays of arrays having tuples (a,b), making sure that the ones in List2
        aren't in List 1
    """
    i = 0
    for sequence in List2:
        for position in sequence:
            if position not in List1[i]:
                List1[i].append(position)
        i += 1

def correctPositions(next, val):
    """
        For each row in next, for each tuple, we increase by val
    """
    copy = next
    next = []
    i = 0
    for row in copy:
        next.append([])
        for position in row:
            new_pos = (position[0] + val, position[1] + val)
            if new_pos not in next[i]: next[i].append(new_pos)
        i += 1
    return next

def getPositionsInsertion(motif, rna):
    """
    Will return an array of array tuples of two positions. Each tuple will be of size len(motif).
    Each tuple represent the positions (STARTING AT 1) where the motifs can be inserted in the rna
    """
    found = re.search(motif[0], rna)
    if found == None: #if we can't match the first sequence, we return None
        return None
    else:
        found = found.start()


    position = [[(found + 1, found + len(motif[0]))]]
    if(len(motif) == 1): #if we can match, and there is only one sequence,
        if found + len(motif[0]) < len(rna):
            next = getPositionsInsertion(motif, rna[(found + 1):]) #  we check if there is another position farrer
            if next == None:
                return position
            else:
                next = correctPositions(next, found + 1)
                position[0].extend(next[0])
                return position
        else: #there is not enough space to find another time the motif, we return the position
            return position
    else: #If there is more then one sequence, we must also match the followings
        next = None
        if found + len(motif[0]) + len(motif[1]) - 1 < len(rna) : #we make sure we can found the rest of the motifs
            next = getPositionsInsertion(motif[1:], rna[found + len(motif[0]):])
        #if we can't place the following sequence, we return None
        if next == None: #if we can't match the rest, we return none
            return None
        else:
            next = correctPositions(next, found + len(motif[0])) #if the following are there, we must correct the positions
            position.extend(next) #if yes then we extend our position array of array with the next
            following = getPositionsInsertion(motif, rna[found + 1 :]) #we must check if the motif can be placed elsewhere
            if following == None:
                return position
            else:
                following = correctPositions(following, found + 1)
                mergePositionsLists(position, following)
                return position

def createMotifsDict(rna, pathDesc):
    """ This function takes as input an 'rna' secondary structure and the path to a folder
    that contains motifs .desc files. It will output a dictionary with key the lengths of the motifs
    and value: A dictionary of : the name of motifs
       and value an array of positions (a,b) where each part can be placed  
    A special key: MAX hold as value the biggest number of sequence in a motif that can be inserted
    """
    if not os.path.isdir(pathDesc):
        print 'there is an error with the path to the DESC folder : %s' % pathDesc
        sys.exit()
    dict = {}
    #We go through every .desc file
    for motifDesc in (x for x in os.listdir(pathDesc) if x[-5:] == '.desc'):
        with open(os.path.join(pathDesc, motifDesc), 'r') as f:
            f.readline() #and we just want the information on the second line
            seq = getArrayMotSequence(f.readline()) #Now we retrieve the array with the sequences in an ordered array
            if seq == None: #if there is an error in the sequence, we skip it AND RECORD IT
                with open("logBadDesc.txt" , "a") as f:
                    f.write('%s\n' % motifDesc)
                continue
            positions = getPositionsInsertion(seq, rna)
            if positions == None: #if we can't place the motif, we go to the next one
                continue
            minLengthSequenceThree(positions, len(rna))#need len(rna) for border cases
            tmpDict = {motifDesc[0:-5] : positions}
            key = str(len(positions))
            if dict == {}:
                dict[key] = tmpDict
                dict["max"] = len(positions)
            elif key in dict:
                dict[key].update(tmpDict) #merge should be ok since the key is the name of the motif and they shoud be processed once only
                dict["max"] = dict["max"] if dict["max"] > len(positions) else len(positions)
            else:
                dict[key] = tmpDict
                dict["max"] = dict["max"] if dict["max"] > len(positions) else len(positions)

    if(len(dict) < 2):
        print "No motif could be inserted.\n"
        sys.exit()

    return dict

def restrain_max_nb_components_in_motif_dict(dict, nb):
    tmpdict = {}
    real_max = 1
    for i in range(nb):
        if str(i+1) in dict:
            tmpdict[str(i + 1)] = dict[str(i + 1)].copy()
            real_max = i + 1 
        else:
            tmpdict[str(i+1)] = []
    tmpdict['max'] = real_max
    return tmpdict

def minLengthSequenceThree(position, lenRNA):
    """"
        We adjust the length of sequences to make sure at least of length 3
        input: array of positions
    """
    for sequences in position:
        tmpSeq = sequences[:]
        for place in tmpSeq: #Use tmpSeq to loop since we remove element from sequences, which make the loop iterator miss items
            if place[1] - place[0] < 2:
                tmp = place
                sequences.remove(place)
                if tmp[0] - tmp[1] == 0: #length 1
                    #We must make sure we don't go out of boundaries nor have duplicates
                    if tmp[0] == 1:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                    elif tmp[0] == 2:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                        if (tmp[0] - 1, tmp[0] + 1) not in sequences: sequences.append((tmp[0] - 1, tmp[0] + 1))
                    elif tmp[0] - lenRNA == 0:
                        if (tmp[0] - 2, tmp[0]) not in sequences: sequences.append((tmp[0] - 2, tmp[0]))
                    elif tmp[0] - lenRNA == -1:
                        if (tmp[0] - 2, tmp[0]) not in sequences: sequences.append((tmp[0] - 2, tmp[0]))
                        if (tmp[0] - 1, tmp[0] + 1) not in sequences: sequences.append((tmp[0] - 1, tmp[0] + 1))
                    else:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                        if (tmp[0] - 1, tmp[0] + 1) not in sequences: sequences.append((tmp[0] - 1, tmp[0] + 1))
                        if (tmp[0] - 2, tmp[0]) not in sequences: sequences.append((tmp[0] - 2, tmp[0]))

                else: #since we consider lengths smaller then 3, it must be of length 2 (no null by construction)
                    if tmp[0] == 1:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                    elif tmp[0] == 2:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                        if (tmp[0] - 1, tmp[0] + 1) not in sequences: sequences.append((tmp[0] - 1, tmp[0] + 1))
                    elif tmp[1] - lenRNA == 0:
                        if (tmp[1] - 2, tmp[1]) not in sequences: sequences.append((tmp[1] - 2, tmp[1]))
                    elif tmp[1] - lenRNA == -1:
                        if (tmp[1] - 2, tmp[1]) not in sequences: sequences.append((tmp[1] - 2, tmp[1]))
                        if (tmp[1] - 1, tmp[1] + 1) not in sequences: sequences.append((tmp[1] - 1, tmp[1] + 1))
                    else:
                        if (tmp[0], tmp[0] + 2) not in sequences: sequences.append((tmp[0], tmp[0] + 2))
                        if (tmp[0] - 1, tmp[0] + 1) not in sequences: sequences.append((tmp[0] - 1, tmp[0] + 1))

def create_motifs_names_dict(motifs_dict):
    motifs_name_dict = {}
    for i in range(motifs_dict['max']):
        motifs_name_dict['MOTIFS%d' % (i+1)] = [x for x in motifs_dict[str(i + 1)]]
    return motifs_name_dict

def create_components_dict(motifs_dict, motifs_names_dict):
    cpts_dict = {}
    for i in range(motifs_dict['max']):
        cpts_dict['CPTS%d' % (i+1)] = []
        for j in range(i + 1):
            cpts_dict['CPTS%dOf%d' % (j+1, i+1)] = []
            for mot in motifs_names_dict['MOTIFS%d' % (i+1)]:
                cpts_dict['CPTS%dOf%d' % (j+1, i+1)].extend([(mot,k,l) for (k,l) in motifs_dict[str(i+1)][mot][j]])
            cpts_dict['CPTS%d' % (j+1)].extend(cpts_dict['CPTS%dOf%d' % (j+1, i+1)])
    return cpts_dict

def add_components_variables_gurobi(m, cpts_dict, vars_dict, max_cpts):
    for i in range(max_cpts):
        for (x,k,l) in cpts_dict['CPTS%d' % (i+1)]:
            cpt_name = "C-%s-%d-%d-%d" % (x,k,l,i+1)
            vars_dict[cpt_name] = m.addVar(vtype = GRB.BINARY, name = cpt_name)
    return None

def add_bases_variables_gurobi(m,BASES, vars_dict):
    for(i, j) in BASES:
        D = "D-%d-%d" % (i, j)
        vars_dict[D] = m.addVar(vtype = GRB.BINARY, name = D)
    return None

def make_cpts_weights_dict(cpts_dict, max_cpts):
    weights_dict = {} 
    for i in range(max_cpts):
        for (x,k,l) in cpts_dict['CPTS%d' % (i+1)]:
            if x not in weights_dict:
                weights_dict[x] = ( l-k + 1 , i + 1)
            elif weights_dict[x][1] < i + 1:
                weights_dict[x] = (weights_dict[x][0] + l - k + 1, i + 1)
    return weights_dict

def constraint_objective(m, cpts_dict, BASES, vars_dict, max_cpts):
    """
            minimize min:N-2*sum{(i,j) in BASES}(1-D[i,j]) - (sum{(x,k,l,u) in SEQE}((l-k+1)*M[x,k,l,u]));
    """
    weights = []
    cpts = []
    weights_dict = make_cpts_weights_dict(cpts_dict, max_cpts)
    for (x,k,l) in cpts_dict['CPTS1']:
        cpts.append(vars_dict["C-%s-%d-%d-1" % (x,k,l)])
        weights.append(-(weights_dict[x][0] ** 2))

    for (u, v) in BASES:
        cpts.append(vars_dict["D-%d-%d" % (u, v)])
        weights.append(10.0)
    lin_expr = LinExpr(weights, cpts)
    m.setObjective(lin_expr, 1) #1 for minimization, -1 for max
    return None

def constraint_extremities_only_overlap(m, cpts_dict, vars_dict, max_cpts, BASES, rna_len):
    for z in range(rna_len):
        arround = []
        over = []
        for i in range(max_cpts):
            arround += [vars_dict["C-%s-%d-%d-%d" % (x,k,l,i+1)] for (x,k,l) in cpts_dict['CPTS%d' % (i+1)] if k < z+1 < l]
            over += [vars_dict["C-%s-%d-%d-%d" % (x,k,l,i+1)] for (x,k,l) in cpts_dict['CPTS%d' % (i+1)] if k == z+1  or z+1 == l]
        bp_arround = [vars_dict["D-%d-%d" % (u,v)] for (u,v) in BASES if (u == z+1) or (v == z+1)]
        if len(arround) + len(bp_arround) > 0:
            weights = [1.0 for y in range(len(arround))] + [0.75 for y in range(len(over))] + [-0.25 for y in range(len(bp_arround))]
            lin_expr = LinExpr(weights, arround + over + bp_arround)
            lin_expr.addConstant(float(len(bp_arround))*0.25)
            m.addConstr(1.0001, GRB.GREATER_EQUAL, lin_expr, "only_one_%d" % (z+1))
    return None

def constraint_component_followed(m, cpts_dict, vars_dict, max_cpts):
    for i in range(1, max_cpts):
        for j in range(i):
            for (x,k,l) in cpts_dict["CPTS%dOf%d" % (j+1, i+1)]:
                next_cpts = [vars_dict["C-%s-%d-%d-%d" % (x2,k2,l2,j+2)] for (x2,k2,l2) in cpts_dict['CPTS%dOf%d' % (j+2, i+1)] if
                        x == x2 and k2 > l + 5]
                m.addConstr(vars_dict["C-%s-%d-%d-%d" % (x,k,l,j+1)], GRB.LESS_EQUAL, LinExpr([1.0 for z in range(len(next_cpts))], next_cpts),
                    "Followed%dOf%d-%s-%d-%d" % (j+1, i+1, x, k, l))
    return None

def constraint_component_preceded(m, cpts_dict, vars_dict, max_cpts):
    for i in range(1, max_cpts):
        for j in range(i):
            for (x,k,l) in cpts_dict["CPTS%dOf%d" % (j+2, i+1)]:
                prev_cpts = [vars_dict["C-%s-%d-%d-%d" % (x2,k2,l2,j+1)] for (x2,k2,l2) in cpts_dict['CPTS%dOf%d' % (j+1, i+1)] if
                        x == x2 and l2 < k - 5]
                m.addConstr(vars_dict["C-%s-%d-%d-%d" % (x,k,l,j+2)], GRB.LESS_EQUAL, LinExpr([1.0 for z in range(len(prev_cpts))], prev_cpts),
                    "Preceded%dOf%d-%s-%d-%d" % (j+2, i+1, x, k, l))
    return None

def constraint_motifs_entirely_inserted(m, cpts_dict,motifs_names_dict, vars_dict, max_cpts):
    for i in range(1, max_cpts):
        for mot in motifs_names_dict['MOTIFS%d' % (i+1)]:
            first = [vars_dict["C-%s-%d-%d-%d" % (x,k,l,1)] for (x,k,l) in cpts_dict["CPTS%dOf%d" % (1, i+1)] if x == mot]
            next_cpts = []
            for j in range(1,i+1):
                next_cpts.extend([vars_dict["C-%s-%d-%d-%d" % (x,k,l,j+1)] for (x,k,l) in cpts_dict["CPTS%dOf%d" % (j+1, i+1)] if x == mot])
            weights = [ i for z in range(len(first))] + [-1.0 for z in range(len(next_cpts))]
            m.addConstr(LinExpr(weights, first + next_cpts), GRB.EQUAL, 0, "AllSeq-%s" % mot)
    return None

def constraint_max_basepairs_removal(m, BASES, max_bp_removal, vars_dict):
    bps = [vars_dict["D-%d-%d" % (i, j)] for (i, j) in BASES]
    m.addConstr(LinExpr([1.0 for z in range(len(bps))], bps), GRB.LESS_EQUAL, len(BASES) * max_bp_removal, "max_nb_bp_removed")
    return None

def constraint_components_surrounded(m, cpts_dict, motifs_names_dict, vars_dict, BASES, max_cpts):
    for i in range(max_cpts):
        for (x,k,l) in cpts_dict["CPTS%d" % (i+1)]:
            if x in motifs_names_dict['MOTIFS1']:
                bps_arround = [vars_dict["D-%d-%d" % (u,v)] for (u,v) in BASES if  k >= u >= k-1 and l <= v <= l + 1 ]
                cpts_arround = [vars_dict["C-%s-%d-%d-%d" % (x2,k2,l2,1)] for (x2,k2,l2) in cpts_dict['CPTS1Of2'] if   l2 == k - 1]
                cpts_arround += [vars_dict["C-%s-%d-%d-%d" % (x2,k2,l2,2)] for (x2,k2,l2) in cpts_dict['CPTS2Of2'] if  k2 == l + 1]
                if len(bps_arround) + len(cpts_arround) > 0:
                    weights = [-1.0 for z in range(len(bps_arround))] + [1.0 for z in range(len(cpts_arround))]
                    lin_expr = LinExpr(weights, bps_arround + cpts_arround)
                    lin_expr.addConstant(len(bps_arround))
                    m.addConstr(lin_expr, GRB.GREATER_EQUAL, vars_dict["C-%s-%d-%d-1" % (x,k,l)], "cpt_surrounded-%s-%d-%d-1" %(x,k,l))
                else:
                    m.addConstr(0, GRB.GREATER_EQUAL, vars_dict["C-%s-%d-%d-1" % (x,k,l)], "cpt_surrounded-%s-%d-%d-1" %(x,k,l))
            else:
                bps_arround = [vars_dict["D-%d-%d" % (u,v)] for (u,v) in BASES if k >= u >= k-1 or l <= u <= l+1 or k >= v >= k-1 or l <= v <= l + 1 ]
                if len(bps_arround) > 0:
                    weights = [-1.0 for z in range(len(bps_arround))]
                    lin_expr = LinExpr(weights, bps_arround)
                    lin_expr.addConstant(len(bps_arround))
                    m.addConstr(lin_expr, GRB.GREATER_EQUAL, vars_dict["C-%s-%d-%d-%d" % (x,k,l, i+1)], "cpt_surrounded-%s-%d-%d-%d" %(x,k,l,i+1))
                else:
                    m.addConstr(0, GRB.GREATER_EQUAL, vars_dict["C-%s-%d-%d-%d" % (x,k,l, i+1)], "cpt_surrounded-%s-%d-%d-%d" %(x,k,l,i+1))
    return None

def constraint_most_one_motif_more_3_cpts(m, cpts_dict, vars_dict, max_cpts):
    motifs_list = []
    for i in range(2,max_cpts):
        motifs_list.extend([vars_dict["C-%s-%d-%d-1" % (x,k,l)] for (x,k,l) in cpts_dict["CPTS1Of%d" % (i+1)]])
        m.addConstr(LinExpr([1.0 for z in range(len(motifs_list))], motifs_list), GRB.LESS_EQUAL, 1, "1MotOf3Seqs")
    return None

def constraint_interior_loops_arround_well_balanced(m, BASES, cpts_dict, vars_dict, motifs_names_dict, rna_len):
    #We want to have, for motifs with two sequences x_1, x_2, for each BP k,l: -n*D_k,k <= (sum x_1 <k - sum x_2 <k -[sum x_2 >l - sum x_1 >l] <= n*D_k,l
    #so if there is a BP, we must have exactly as much first seq before not filled as second seq after not filled 
    for (u,v) in BASES:
        for mot in motifs_names_dict['MOTIFS2']:
            first = [vars_dict["C-%s-%d-%d-1" % (x,k,l)] for (x,k,l) in cpts_dict['CPTS1Of2'] if mot == x and (l < u or k > v)]
            second = [vars_dict["C-%s-%d-%d-2" % (x,k,l)] for (x,k,l) in cpts_dict['CPTS2Of2'] if mot == x and (l < u or k > v)]
            weights = [1.0 for z in range(len(first))] + [-1.0 for z in range(len(second))]
            m.addConstr(LinExpr([-rna_len], [vars_dict["D-%d-%d" % (u,v)]]), GRB.LESS_EQUAL, LinExpr(weights, first + second), "2_cpts_consistend_left-%d-%d" % (u,v))
            m.addConstr(LinExpr([rna_len], [vars_dict["D-%d-%d" % (u,v)]]), GRB.GREATER_EQUAL, LinExpr(weights, first + second), "2_cpts_consistend_left-%d-%d" % (u,v))
    return None

def constraint_3_way(m, cpts_dict, vars_dict, BASES, motifs_names_dict, rna_len):
    #Make sure all compts are close by the same "branch"
    for (k,l) in BASES:
        #Three way junctions constraints, must "inside" the same BPs
        first_part =  [vars_dict["C-%s-%d-%d-1" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS1Of3'] if i >= k and j <= l] 
        second_part = [vars_dict["C-%s-%d-%d-2" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS2Of3'] if i >= k and j <= l] 
        third_part =  [vars_dict["C-%s-%d-%d-3" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS3Of3'] if i >= k and j <= l] 
        coeffs = [2.0 for x in range(len(first_part))]
        coeffs.extend([-1.0 for x in range(len(second_part) + len(third_part))])
        tot_part = []
        tot_part.extend(first_part)
        tot_part.extend(second_part)
        tot_part.extend(third_part)
        m.addConstr(LinExpr(coeffs, tot_part), GRB.LESS_EQUAL, LinExpr( [rna_len] ,[vars_dict["D-%d-%d" % (k,l)]]), "+3_way_between_D-%d-%d" % (k,l))
        m.addConstr(LinExpr(coeffs, tot_part), GRB.GREATER_EQUAL, LinExpr( [-rna_len] ,[vars_dict["D-%d-%d" % (k,l)]]), "-3_way_between_D-%d-%d" % (k,l))

def constraint_4_way(m, cpts_dict, vars_dict, BASES, motifs_names_dict, rna_len):
    for (k,l) in BASES:
        first_part =  [vars_dict["C-%s-%d-%d-1" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS1Of4'] if i >= k and j <= l] 
        second_part = [vars_dict["C-%s-%d-%d-2" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS2Of4'] if i >= k and j <= l] 
        third_part =  [vars_dict["C-%s-%d-%d-3" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS3Of4'] if i >= k and j <= l] 
        fourth_part = [vars_dict["C-%s-%d-%d-4" % (x,i,j)] for (x,i,j) in cpts_dict['CPTS4Of4'] if i >= k and j <= l] 
        coeffs = [3.0 for x in range(len(first_part))]
        coeffs.extend([-1.0 for x in range(len(second_part) + len(third_part) + len(fourth_part))])
        tot_part = []
        tot_part.extend(first_part)
        tot_part.extend(second_part)
        tot_part.extend(third_part)
        tot_part.extend(fourth_part)
        m.addConstr(LinExpr(coeffs, tot_part), GRB.LESS_EQUAL, LinExpr( [rna_len] ,[vars_dict["D-%d-%d" % (k,l)]]), "+4_way_between_D-%d-%d" % (k,l))
        m.addConstr(LinExpr(coeffs, tot_part), GRB.GREATER_EQUAL, LinExpr( [-rna_len] ,[vars_dict["D-%d-%d" % (k,l)]]), "-4_way_between_D-%d-%d" % (k,l))

def constraint_no_lonely_bp(m, vars_dict, BASES, rna_len):
    #First border cases
    #i == 1
    i = 1
    before_after = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i + 1 or v == i + 1]
    between = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i or v == i]
    weights = [-1.0 for z in range(len(before_after))] + [1.0 for z in range(len(between))]
    lin_expr = LinExpr(weights, before_after + between)
    if len(between) == 1:
        lin_expr.addConstant(len(before_after))
    else:
        lin_expr.addConstant(len(before_after) + 1)
    m.addConstr(1, GRB.LESS_EQUAL, lin_expr, "stack_bp_%d" % i)
    #i == rna_len 
    i = rna_len 
    before_after = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i - 1 or v == i - 1]
    between = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i or v == i]
    weights = [-1.0 for z in range(len(before_after))] + [1.0 for z in range(len(between))]
    lin_expr = LinExpr(weights, before_after + between)
    if len(between) == 1:
        lin_expr.addConstant(len(before_after))
    else:
        lin_expr.addConstant(len(before_after) + 1)
    
    for i in range(2, rna_len):
        before_after = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i-1 or u == i + 1 or v == i-1 or v == i + 1]
        between = [vars_dict['D-%d-%d' % (u,v)] for (u,v) in BASES if u == i or v == i]
        weights = [-1.0 for z in range(len(before_after))] + [1.0 for z in range(len(between))]
        lin_expr = LinExpr(weights, before_after + between)
        if len(between) == 1:
            lin_expr.addConstant(len(before_after))
        else:
            lin_expr.addConstant(len(before_after) + 1)
        m.addConstr(1, GRB.LESS_EQUAL, lin_expr, "stack_bp_%d" % i)

def constraint_cover_more_bp(model, cpts_dict, vars_dict, BASES, motifs_names_dict ):
    for mot in motifs_names_dict['MOTIFS2']:
        for (x,k,l) in ( (x2,k2,l2) for (x2,k2,l2) in cpts_dict['CPTS1Of2'] if x2 == mot):
            for (y,m,n) in ( (x2,k2,l2) for (x2,k2,l2) in cpts_dict['CPTS2Of2'] if x2 == mot and k2 > l + 3):
                coverage = list(range(k,l+1)) + list(range(m,n+1))
                bases_inside = [(u,v) for (u,v) in BASES if u in coverage and v in coverage]
                bases = [(u,v) for (u,v) in BASES if u in coverage or v in coverage and (u,v) not in bases_inside]
                if  len(coverage)-1  <= 2*len(bases_inside)+len(bases):
                    model.addConstr(LinExpr([1,1], [vars_dict['C-%s-%d-%d-1' % (mot, k,l)], vars_dict['C-%s-%d-%d-2' % (mot, m, n)]]), GRB.LESS_EQUAL, 1)

def gurobi_create_model(motifs_dict, secStructPos, rna, max_bp_removal):
    """
        We create the sets that will be used in our model.
        return m, [vars_dict, cpts_dict, motifs_names_dict, BASES]
    """
    rna_len = len(rna)

    max_cpts = motifs_dict['max']

    motifs_names_dict = create_motifs_names_dict(motifs_dict)

    cpts_dict = create_components_dict(motifs_dict, motifs_names_dict)

    BASES = [(u, v) for (u, v) in secStructPos]

    # WE CREATE OUR MODEL
    m = Model("toMinimize")
    
    #We need a dict to keep track of the variables "names"
    vars_dict = {}

    #Create Variables
    add_components_variables_gurobi(m, cpts_dict, vars_dict, max_cpts)

    add_bases_variables_gurobi(m,BASES, vars_dict)

    m.update()

    constraint_objective(m, cpts_dict, BASES, vars_dict, max_cpts)
    constraint_cover_more_bp(m, cpts_dict, vars_dict, BASES, motifs_names_dict)

    """===========================
            CONSTRAINTS
    ==========================="""
    constraint_no_lonely_bp(m, vars_dict, BASES, rna_len)
    constraint_extremities_only_overlap(m, cpts_dict, vars_dict, max_cpts,BASES, rna_len)

    constraint_component_followed(m, cpts_dict, vars_dict, max_cpts)
    constraint_component_preceded(m, cpts_dict, vars_dict, max_cpts)

    constraint_motifs_entirely_inserted(m, cpts_dict,motifs_names_dict, vars_dict, max_cpts)

    constraint_max_basepairs_removal(m, BASES, max_bp_removal, vars_dict)

    constraint_components_surrounded(m, cpts_dict, motifs_names_dict, vars_dict, BASES, max_cpts)

    constraint_most_one_motif_more_3_cpts(m, cpts_dict, vars_dict, max_cpts)

    constraint_interior_loops_arround_well_balanced(m, BASES, cpts_dict, vars_dict, motifs_names_dict, rna_len)

    constraint_3_way(m, cpts_dict, vars_dict, BASES, motifs_names_dict, rna_len)

    if max_cpts > 3:
        constraint_4_way(m, cpts_dict, vars_dict, BASES, motifs_names_dict, rna_len)

    m.update()
    return m, [vars_dict, cpts_dict, motifs_names_dict, BASES]

def gurobi_find_all_optimal_solutions(motifs_dict, sec_struct_pos, rna, max_bp_removal,max_nb_sols):
    """Solve for all optimal solutions given a motifs_dict, sec_struct_pos, rna, max_bp_removal
        will output a tuple with elements:
            1: a list of of list if the different "C-%s-%d-%d-%d" == 1 in each solutions
            2: the gurobi model after all solutions where generated (see item 4)
            3: the list model_dict generate for the given problem gurobi_create_model
            4: the list of lin_expr added to the model to generate all optimal solutions
    """
    nb_sols = 1 
    list_sols = []
    list_of_lin_expr = []

    m, model_list_of_dicts = gurobi_create_model(motifs_dict, sec_struct_pos, rna, max_bp_removal)
    m.setParam("Threads",NUMBER_THREADS)
    m.setParam("MIPFocus",3)
    m.setParam("MIPGap", 0)
    m.setParam("Method", 3)
    m.update()
    m.optimize()

    if m.Status != GRB.Status.OPTIMAL:
            print '-------------------------'
            print ' NO SOLUTIONS!'
            sys.exit(0)

    min_val = m.ObjVal
    while True:  #we loop for all optimal solutions
        nb_sols += 1
        in_sol = []
        not_in_sol = []
        list_sols.append([])

        for x in m.getVars():
            if round(x.X)  == 1:
                in_sol.append(x)
                list_sols[-1].append(x.VarName)
            else:
                not_in_sol.append(x)

        if nb_sols > max_nb_sols:
            break

        weights = [-1.0 for z in range(len(in_sol))]
        weights.extend([1.0 for z in range(len(not_in_sol))])
        tot_sols = in_sol + not_in_sol
        lin_expr = LinExpr(weights, tot_sols)
        list_of_lin_expr.append(lin_expr)
        m.addConstr(lin_expr, GRB.GREATER_EQUAL, -len(in_sol) + 1)
        m.update()
        print '\n#########################\n'
        m.optimize()
        if not m.Objval or m.ObjVal > min_val + 0.5:
            print "\nAll Optimal Solutions have been found\nNext best score: %s\n" % m.Objval 
            break

    return m, model_list_of_dicts, list_sols, list_of_lin_expr, min_val

def validate_rna_seq(rna_seq):
    rna_seq = rna_seq.strip().upper()
    if not all(x in 'ACGU' for x in rna_seq):
        help(rna_seq)
        sys.exit(1)
    return rna_seq

def validate_sec_struct(sec_struct_info):
    l_struct_positions = []

    if os.path.isfile(sec_struct_info):
        l_sec_struct = [x.strip() for x in open(sec_struct_info)
                        if x]
    else:
        l_sec_struct = [sec_struct_info.strip()]

    for sec_struct in l_sec_struct:
        sec_struct_positions = []
        if not all(x in '(.)' for x in sec_struct):
            help(sec_struct=1)
            sys.exit(1)
        left = []
        for i,pos in enumerate(sec_struct):
            if pos == '(':
                left.append(i)
            elif pos == ')':
                try:
                    l = left.pop()
                except IndexError:
                    help(sec_struct=1)
                    sys.exit(1)
                sec_struct_positions.append(
                    (l+1,i+1))
        if left:
            help(sec_struct=1)
            sys.exit(1)
        l_struct_positions.append((sec_struct,sec_struct_positions))
    return l_struct_positions

def validate_max_bp_removal(max_bp_removal):
    max_bp_removal = float(max_bp_removal.strip())
    if not 0 <= max_bp_removal <= 1.0:
        help(max_bp_removal=1)
        sys.exit(1)
        
    return max_bp_removal

def validate_max_components(max_components):
    max_components = int(max_components.strip())
    if max_components < 3:
        help(max_components=1)
        sys.exit(1)
    return max_components

def validate_max_nb_sols(max_nb_sols):
    max_nb_sols = int(max_nb_sols.strip())
    if max_nb_sols < 1:
        help(max_nb_sols=1)
        sys.exit(1)
    return max_nb_sols

def remove_deleted(sec_struct,list_sols):
    ss = list(sec_struct)
    for x in list_sols:
        x = x.split('-')
        if x[0] == 'D':
            ss[int(x[1])-1] = '.'
            ss[int(x[2])-1] = '.'
    return ''.join(ss)

def help(rna_seq='',
         sec_struct='',
         max_bp_removal='',
         max_components='',
         required_arg='',
         max_nb_sols=''):
    if rna_seq:
        print 'The RNA sequence should only contain the characters "ACGU"\n'
    if sec_struct:
        print """The secondary structure should only contain the characters "(.)"
            and be a well balanced parenthese equation\n"""
    if max_bp_removal:
        print "The max number of base pairs should be between 0 and 1\n"
    if max_components:
        print "The number of components should be greater than 3\n"
    if required_arg:
        print "A required argument is missing\n"
    if max_nb_sols:
        print "The maximal number of solutions should be bigger than 0\n"

    print """ RNAMoIP v1.1
        ARGUMENTS 
          REQUIRED:
            -s The rna Sequence 
            -ss The rna Secondary Structure (no pseudoknots) 
            -d the path to the .Desc files
          OPTIONAL
            -r (default 0.3) 
               the percentage of basepairs that can be Removed 
            -c (default 4)
                the max nb of Components in motifs 
            -m_sols (default 1)
                maximal number of optimal solutions to output
        e.g.
            gurobi.sh RNAMoIP.py 'GGGCGGCCUUCGGGCUAGACGGUGGGAGAGGCUUCGGCUGGUCCACCCGUGACGCUC' '((((((((....))))..((((..(((..(((....)))..)))..))))...))))' 'No_Redondance_DESC' 0.3 4 > my_output.txt
    """

if __name__ == '__main__':
    rna_seq = ''
    l_sec_struct_positions = []
    path_desc = ''
    max_nb_sols = 1
    max_bp_removal = 0.3
    max_components = 4

    opts = sys.argv
    for i,option in enumerate(opts):
        if option == '-s':
            rna_seq = validate_rna_seq(opts[i+1])
        elif option == '-ss':
            l_sec_struct_positions = validate_sec_struct(opts[i+1])
        elif option == '-d':
            path_desc = opts[i+1]
        elif option == '-r':
            max_bp_removal = validate_max_bp_removal(opts[i+1])
        elif option == '-c':
            max_components = validate_max_components(opts[i+1])
        elif option == '-m_sols':
            max_nb_sols = validate_max_nb_sols(opts[i+1])
        elif option in ('-h','-help'):
            print help()
            sys.exit(0)
    if not all((rna_seq,
                l_sec_struct_positions,
                path_desc)):
        help(required_arg=1)
        sys.exit(1)

    """=============================================================================================
    Generate models and sols
    ============================================================================================="""
    motifs_dict = createMotifsDict(rna_seq, path_desc)
    motifs_dict = restrain_max_nb_components_in_motif_dict(motifs_dict, max_components)

    all_sols = []
    for sec_struct,sec_struct_positions in l_sec_struct_positions:
        m, model_list_of_dicts, list_sols, list_of_lin_expr, min_val =  gurobi_find_all_optimal_solutions(motifs_dict,
                                                                                             sec_struct_positions,
                                                                                             rna_seq,
                                                                                             max_bp_removal,
                                                                                             max_nb_sols)
        all_sols.append((list_sols, min_val,sec_struct))
    all_sols = sorted(all_sols,key=lambda x:x[1])
    min_val = all_sols[0][1]

    for list_sols,val,sec_struct in all_sols:
        if val > min_val:
            break
        print "Solution for the secondary structure:"
        print "\t%s" %  sec_struct
        for i in range(len(list_sols)):
            print
            print 'Optimal solution nb: ',i + 1
            print 'Corrected secondary structure:'
            print "\t%s" % remove_deleted(sec_struct,list_sols[i])
            for name in list_sols[i]:
                print '\t', name

    print "\nThe optimal solutions has as value:\n\t%s" % min_val,
