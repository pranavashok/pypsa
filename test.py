from psa.LearnPSA import LearnPSA

L = LearnPSA(10, 10, 3, ["c", "d", "e"])
L.learn_sample("cdecdcdecd")
L.learn_sample("dedecdcdec")
#print((1+L.e2)*L.gamma_min, 1+3*L.e2, (1-L.e1)*L.e0)
L.print_tree()
states, transition = L.generate_psa()
L.generate_run(states, transition, 10)
#print(L._P2('d', 'c'))