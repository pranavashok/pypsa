from psa.LearnPSA import LearnPSA

L = LearnPSA(0.8, 50, 9, ["a", "b", "c", "d", "e", "f", "g"])
L.learn_sample("cdecdcdecd")
L.learn_sample("dedecdcdec")
L.print_tree()