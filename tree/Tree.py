class Tree(object):
    def __init__(self, data):
        self.parent = None
        self.children = []
        self.data = data
    
    def insert(self, n):
        n.parent = self
        self.children.append(n)

    def bfs(self, query):
        bfsqueue = []
        for c in self.children:
            bfsqueue.append(c)
        while len(bfsqueue) != 0:
            for e in bfsqueue:
                if e.data == query:
                    return e
                else:
                    for c in e.children:
                        if c.data == query:
                            return c
                        else:
                            bfsqueue.append(c)

    def path_to_root(self):
        node = self
        path = []
        while node:
            path.append(node.data)
            node = node.parent

        path.reverse()

        for i in path:
            print(i)