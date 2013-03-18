from tree.Tree import Tree
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

    def learn_sample(self, s):
        self.sample.append(s)
        self.PST = self._learn()

    def _P1(self, s):
        '''
        P(s) is roughly the relative number of times s appears in the sample
        This implementation of P is slightly modified. It divides with |r| - L 
        for each string r in sample and each |r| does not need to be equal to l
        '''
        p = 0
        for r in self.sample:
            p += r[self.L-len(s)+1:len(r)-1].count(s)
        p = p/(len(self.sample)*(len(r)-self.L))
        return p

    def _P2(self, sigma, s):
        '''
        P(sigma|s) is roughly the relative number of times sigma appears after s
        '''
        countssigma = 0
        counts = 0
        for r in self.sample:
            countssigma += r[self.L-len(s+sigma)+1:len(r)-1].count(s+sigma)
            counts += r[self.L-len(s)+1:len(r)-1].count(s)
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
                    if tree.data[0] == '0':
                        missing_children.append(sigma)
                    else:
                        missing_children.append(sigma+tree.data[0])
            else:
                return tree

            for child in tree.children:
                child = self._add_missing_children(child)
                if child.data[1] == 1:
                    missing_children.remove(child.data[0])

            for s in missing_children:
                tree = tree.insert(Tree([s, 0]))

            return tree

    def _learn(self):
        T = Tree(["0", 1])
        S = []
        removed_from_S = []
        for sigma in self.Sigma:
            if self._P1(sigma) >= (1-self.e1)*self.e0:
                S.append(sigma)

        while len(S) > 0:
            s = S.pop()

            for sigma in self.Sigma:
                '''suffix(s) = s[1:]'''
                if s[1:] == "":
                    if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min):
                        T = T.insert(Tree([s, 1]))
                        break
                else:
                    if self._P2(sigma, s[1:]) != 0:
                        if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min) and ((self._P2(sigma, s)/self._P2(sigma, s[1:])) > 1+3*self.e2):
                            '''
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
                                if removed_from_S.count(s[i:]) == 0:
                                    S.append(s[i:]) #Insert only uniques
                                i = i - 1
                            break
                        else:
                            if self._P2(sigma, s) >= (1+self.e2)*self.gamma_min:
                                '''
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
                    if self._P1(sigma+s) >= (1-self.e1)*self.e0:
                        S.append(sigma+s)

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
            print(e.data[0])
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
        for state in psa:
            for sigma in self.Sigma:
                for i in range(0, len(state+sigma)-1):
                    #If state+sigma or it's suffix is present
                    if psa.count((state+sigma)[i:]) == 1:
                        transition[(state, sigma)] = ((state+sigma)[i:], self._P2(sigma, state))
                        flag = 0
                        break
                    else:
                        flag = 1
                if flag == 1:
                    transition[(state, sigma)] = (state+sigma, 0)

        return psa, transition

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

    def generate_run(self, states, transition, N):
        run = ""
        first = self._first_transition()
        run += first
        cur_state = first
        while len(run) <= N:
            print(cur_state)
            po = []
            T = {}
            for sigma in self.Sigma:
                po.append(transition[(cur_state, sigma)][1])
            prob = numpy.array(po)
            cumprob = numpy.cumsum(prob)
            i = 0
            for sigma in self.Sigma:
                T[sigma] = cumprob[i]
                i += 1
            print(T)
            r = random.randrange(1, 999)/1000
            for sigma in self.Sigma:
                if T[sigma]-r >= 0 and T[sigma] >= 0:
                    #If there exists some non-zero transition from sigma|cur_state
                    for s in self.Sigma:
                        flag = 0
                        if transition[(transition[(cur_state, sigma)][0], s)][1] != 0:
                            flag = 1
                            break
                    if flag == 1:
                        run += sigma
                        cur_state = transition[(cur_state, sigma)][0]
                        break
        return run