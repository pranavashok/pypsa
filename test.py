from tree.Tree import Tree

root = Tree(1)

n = Tree(2)

x = Tree(3)
n.insert(x)

x = Tree(4)
n.insert(x)

x = Tree(5)
root.insert(x)

root.insert(n)

result = root.bfs(4)

result.path_to_root()