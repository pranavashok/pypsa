from tree.Tree import Tree
import math

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
            p += r.count(s)/(len(r)-self.L)
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

        if counts > 0:
            p = countssigma/counts
        else:
            p = 0
        return p

    def _add_missing_children(self, tree):
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
            tree = tree.insert(Tree((s, 0)))

        return tree

    def _learn(self):
        T = Tree(("0", 1))
        S = []
        Sigma = self.Sigma
        for sigma in self.Sigma:
            if self._P1(sigma) >= (1-self.e1)*self.e0:
                S.append(sigma)

        while len(S) > 0:
            s = S.pop()

            for sigma in Sigma:
                '''suffix(s) = s[1:]'''
                if s[1:] == "":
                    if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min):
                        T = T.insert(Tree((s, 1)))
                        break
                else:
                    if (self._P2(sigma, s) >= (1+self.e2)*self.gamma_min) and ((self._P2(sigma, s)/self._P2(sigma, s[1:])) > 1+3*self.e2):
                        '''
                        This will insert to the next available suffix. Can it be done this way?
                        '''
                        for i in range(0, self.L):
                            x = T.bfs(s[i:])
                            if x:
                                x = x.insert(Tree((s, 1)))
                                #T = T.bfs(s[i:]).insert(Tree((s, 1)))
                                break

                        break

            if len(s) < self.L:
                for sigma in self.Sigma:
                    if self._P2(sigma, s) >= (1-self.e1)*self.e0:
                        S.append(sigma+s)

        _T = T
        #_T = self._add_missing_children(_T)
        
        '''print(T.data[0])
        for c in T.children:
            print(c.data[0], c.children)
        bfsqueue = []
        for c in _T.children:
            print(c.children)
            bfsqueue.append(c)
        while len(bfsqueue) != 0:
            for e in bfsqueue:
                for d in e.children:
                    bfsqueue.append(d)
        for e in bfsqueue:
            print(e.data[0])'''

        return _T

    def print_tree(self):
        bfsqueue = []
        for c in self.PST.children:
            bfsqueue.append(c)
        while len(bfsqueue) > 0:
            e = bfsqueue.pop()
            print(e.data[0])
            for c in e.children:
                bfsqueue.append(c)