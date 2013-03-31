from __future__ import division
from ..tree.Tree import Tree
import math
import random
import numpy
import array

class LearnPSA(object):
    def __init__(self, e, n, L, Sigma):
        self.e = e
        self.n = n
        self.L = L
        self.Sigma = Sigma
        self.e2 = e/(48*L)
        self.gamma_min = e/(48*L*len(Sigma))
        self.e0 = e/(2*n*L*math.log(48*L*len(Sigma)/e, 10))
        self.e1 = e*math.log(48*L*len(Sigma)/e, 10)/(9216*L*len(Sigma))
        self.sample = []
        self.PST = None
    
    def _count(self, sample, seq):
        count = 0
        flag = 1
        for i in xrange(0, len(sample)-len(seq)+1):
            j = i
            for element in seq:
                if sample[j] == element:
                    flag = 0
                else:
                    flag = 1
                    break
                j = j + 1
            if flag == 0:
                count += 1
        return count

    def _remove(self, l, subl):
        flag = 1
        for i in xrange(0, len(l)-len(subl)+1):
            j = i
            for element in subl:
                if l[j] == element:
                    flag = 0
                else:
                    flag = 1
                    break
                j = j + 1
            if flag == 0:
                for j in xrange(i, i+len(subl)):
                    l.remove(l[i])

    def learn_sample(self, s):
        split = s.split(" ")
        self.sample.append(split)
        self.PST = self._learn()

    def _P1(self, s):
        u'''
        P(s) is roughly the relative number of times s appears in the sample
        This implementation of P is slightly modified. It divides with |r| - L 
        for each string r in sample and each |r| does not need to be equal to l
        '''
        p = 0
        for r in self.sample:
            p += self._count(r[self.L-len(s)+1:len(r)-1], s)
        p = p/(len(self.sample)*(len(r)-self.L))
        return p

    def _P2(self, sigma, s):
        u'''
        P(sigma|s) is roughly the relative number of times sigma appears after s
        '''
        countssigma = 0
        counts = 0
        ssigma = []
        for e in s:
            ssigma.append(e)
        ssigma.append(sigma)
        for r in self.sample:
            countssigma += self._count(r[self.L-len(ssigma)+1:len(r)-1], ssigma)
            counts += self._count(r[self.L-len(s)+1:len(r)-1], s)
        if counts > 0:
            p = countssigma/counts
        else:
            p = 0
        return p

    def _PI(self):
        PI = {}
        for sigma in self.Sigma:
            count = 0
            for r in self.sample:
                if r[0] == sigma:
                    count += 1
            PI[sigma] = count/len(self.sample)
        return PI


    def _add_missing_children(self, tree):
        #Filter out leaves
        if len(tree.children) == 0:
            return tree
        else:
            missing_children = []
            if tree.data[1] == 1:
                for sigma in self.Sigma:
                    if tree.data[0] == u'0':
                        missing_children.append(sigma)
                    else:
                        newnode = [sigma]
                        newnode.append(tree.data[0])
                        missing_children.append(newnode)
            else:
                return tree

            for child in tree.children:
                child = self._add_missing_children(child)
                if child.data[1] == 1:
                    self._remove(missing_children, child.data[0])

            for s in missing_children:
                tree = tree.insert(Tree([s, 0]))

            return tree

    def _learn(self):
        T = Tree([u"0", 1])
        S = []
        removed_from_S = []
        for sigma in self.Sigma:
            if self._P1(sigma) >= (1-self.e1)*self.e0:
                S.append([sigma])

        while len(S) > 0:
            s = S.pop()
            for sigma in self.Sigma:
                u'''suffix(s) = s[1:]'''
                if len(s) == 1:
                    if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min):
                        T = T.insert(Tree([s, 1]))
                        break
                else:
                    if self._P2(sigma, s[1:]) != 0:
                        if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min) and ((self._P2(sigma, s)/self._P2(sigma, s[1:])) > 1+3*self.e2):
                            u'''
                            This will insert all suffixes uptil s
                            '''
                            i = len(s)-1
                            while i >= 0:
                                x = T.bfs(s[i:])
                                if x == None:
                                    parent = T.bfs(s[i+1:])
                                    parent = parent.insert(Tree([s[i:], 1]))
                                i = i - 1
                            i = len(s)-1
                            while i > 0:
                                if s[i:][0] == u" ":
                                    continue
                                if removed_from_S.count(s[i:]) == 0:
                                    S.append(s[i:]) #Insert only uniques
                                i = i - 1
                            break
                        else:
                            if self._P2(sigma, s) >= (1+self.e2)*self.gamma_min:
                                u'''
                                This will insert all suffixes uptil s
                                '''
                                i = len(s)-1
                                while i >= 0:
                                    x = T.bfs(s[i:])
                                    if x == None:
                                        parent = T.bfs(s[i+1:])
                                        parent = parent.insert(Tree([s[i:], 1]))
                                    i = i - 1
                                break

            if len(s) < self.L:
                for sigma in self.Sigma:
                    sigmas = [sigma]
                    sigmas.append(s)
                    if self._P1(sigmas) >= (1-self.e1)*self.e0:
                        S.append(sigmas)

            removed_from_S.append(s)

        _T = T
        _T = self._add_missing_children(_T)
        #_T = self._compute_gamma_s_sigma(_T)

        return _T

    def _compute_gamma_s_sigma(self, tree):
        s = tree.data[0]
        gamma_s_sigma = {}
        for child in tree.children:
            child = self._compute_gamma_s_sigma(child)
        for sigma in self.Sigma:
            gamma_s_sigma[sigma] = self._P2(sigma, (s[1:])[::-1])*(1-len(self.Sigma)*self.gamma_min)+self.gamma_min
            tree.data.append(gamma_s_sigma)
        return tree

    def print_tree(self):
        bfsqueue = []
        for c in self.PST.children:
            bfsqueue.append(c)
        while len(bfsqueue) > 0:
            e = bfsqueue.pop()
            print e.data[0]
            for c in e.children:
                bfsqueue.append(c)

    def _get_pst_states(self):
        states = []
        bfsqueue = []
        for c in self.PST.children:
            bfsqueue.append(c)
        while len(bfsqueue) > 0:
            e = bfsqueue.pop()
            states.append(e.data[0])
            for c in e.children:
                bfsqueue.append(c)
        return states
    
    def generate_psa(self):
        psa = []
        states = self._get_pst_states()
        for state in states:
            psa.append(state[::-1])

        psa.sort()
        transition = {}
        nextstate = {}
        for state in psa:
            for sigma in self.Sigma:
                transition[(state, sigma)] = self._P2(sigma, state)
                if transition[(state, sigma)] > 0:
                    for i in xrange(0, len(state+sigma)-1):
                        #If state+sigma or it's suffix is present
                        if psa.count((state+sigma)[i:]) == 1:
                            nextstate[(state, sigma)] = (state+sigma)[i:]
                            break

        return psa, transition, nextstate

    def _first_transition(self):
        PI = self._PI()
        pi = []
        for sigma in self.Sigma:
            pi.append(PI[sigma])
        prob = numpy.array(pi)
        cumprob = numpy.cumsum(prob)
        i = 0
        for sigma in self.Sigma:
            PI[sigma] = cumprob[i]
            i += 1
        r = random.randrange(0, 1000)/1000
        for sigma in self.Sigma:
            if PI[sigma]-r >= 0 and PI[sigma] >= 0:
                return sigma

    def generate_run(self, states, transition, nextstate, N):
        run = u""
        first = self._first_transition()
        run += first
        cur_state = first
        while len(run) <= 2*N:
            po = []
            T = {}
            for sigma in self.Sigma:
                po.append(transition[(cur_state, sigma)])
            prob = numpy.array(po)
            cumprob = numpy.cumsum(prob)

            i = 0
            for sigma in self.Sigma:
                T[sigma] = cumprob[i]
                i += 1

            r = random.randrange(1, 999)/1000
            for sigma in self.Sigma:
                if T[sigma]-r >= 0 and T[sigma] >= 0:
                    run += u" "+sigma
                    #Find a next possible state which has some non-zero symbol probability
                    next_possible_state = nextstate[(cur_state, sigma)]
                    flag = 0
                    for sigma in self.Sigma:
                        if transition[(next_possible_state, sigma)] > 0:
                            flag = 1
                            break
                    if flag == 1:
                        cur_state = next_possible_state
                        break
                    else:
                        #just setting cur_state to first in case we reach a dead end.
                        cur_state = first
        return run