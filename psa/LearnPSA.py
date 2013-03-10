from ..tree.Tree import Tree
from math import math

class LearnPSA(object):
    def __init__(self, e, n, L, Sigma):
        self.e = e
        self.n = n
        self.L = L
        self.Sigma = Sigma
        self.e2 = e/(48*L)
        self.gamma_min = e2/Sigma
        self.e0 = e/(2*n*L*log(1/gamma_min))
        self.e1 = (e2*gamma_min)/(8*n*e0)

    def _P1(self, s):
        '''
        P(s) is roughly the relative number of times s appears in the sample
        This implementation of P is slightly modified. It divides with |r| - L 
        for each string r in sample and each |r| does not need to be equal to l
        '''
        p = 0
        for r in self.sample:
            p += r.count(s)/(len(s)-L)
        p = p/len(self.sample)

        return p
    
    def _P2(self, sigma, s):
        '''
        P(sigma|s) is roughly the relative number of times sigma appears after s
        '''
        countssigma = 0
        counts = 0
        for r in self.sample:
            countssigma += r.count(s+sigma)
            counts += r.count(s)
        
        p = countssigma/counts
        
        return p

    def learn(self):
        T = Tree(0)
        S = []
        for sigma in Sigma:
            if _P1(sigma) >= (1-e1)*e0:
                S.append(Tree(sigma))

        while len(S) > 0:
            s = S.pop()

            for sigma in Sigma:
                if (_P2(sigma, s.data) >= (1+e2)*gamma_min) && (_P2(sigma, s.data)/_P2(sigma, s.data[1:])) > 1+3*e2): '''suffix(s) = s[1:]'''
                    s = bfs(s.data[1:]).insert(s)
                    '''
                    Is this the correct place in the tree to insert s? What if suffix(s) doesn't exist?
                    '''
                    iter = s.parent
                    while iter:
                        S.append(iter) ''' Should only unique iter be appended? '''
                        iter = iter.parent

            if len(s) < L:
                for sigma in Sigma:
                    if P2(sigma, s.data) >= (1-e1)/e0:
                        S.append(Tree(sigma+s))

        _T = T

