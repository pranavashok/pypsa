from psa.LearnPSA import LearnPSA
from operator import itemgetter

L = LearnPSA(0.2, 10, 3, ["a", "b", "c", "d", "e", "f", "g"])
L.learn_sample("cdecdcdecd")
L.learn_sample("dedecdcdec")
#print((1+L.e2)*L.gamma_min, 1+3*L.e2, (1-L.e1)*L.e0)
#L.print_tree()
states, transition = L.generate_psa()

#Run 100 times and print the frequency of differenct outcomes
run = []
for i in range(0,100):
	run.append(L.generate_run(states, transition, 10))

unique = {}
for item in run:
    unique[item] = unique.get(item, 0) + 1

freq_table = sorted(unique.items(), key=itemgetter(1), reverse=True)
for item in freq_table:
	print(item[0], item[1])