class Tree(object):
    def __init__(self, data):
        self.parent = None
        self.children = []
        self.data = data
    
    def insert(self, n):
        n.parent = self
        self.children.append(n)
        return self

    def bfs(self, query):
        bfsqueue = []
        for c in self.children:
            bfsqueue.append(c)
        while len(bfsqueue) > 0:
            e = bfsqueue.pop()
            if e.data[0] == query:
                return e
            else:
                for c in e.children:
                    if c.data[0] == query:
                        return c
                    else:
                        bfsqueue.append(c)