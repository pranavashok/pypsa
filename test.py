from psa.LearnPSA import LearnPSA

L = LearnPSA(100, 50, 10, ["a", "b", "c", "d", "e", "f", "g"])
L.learn_sample("cdecdcdecda")
L.learn_sample("dedecdcadec")
#print((1+L.e2)*L.gamma_min, 1+3*L.e2, (1-L.e1)*L.e0)
L.print_tree()