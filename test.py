from tree.Tree import Tree

root = Tree(1)

n = Tree(2)

x = Tree(3)
n = n.insert(x)

x = Tree(4)
n = n.insert(x)

x = Tree(5)
root = root.insert(x)

root = root.insert(n)

result = root.bfs(4)

result.path_to_root()