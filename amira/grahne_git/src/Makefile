flags = -Wall -O3 -std=gnu++17

all:: fim_all fim_closed fim_maximal
fim_all: fpgrowth.o fp_tree.o buffer.o fsout.h fp_node.o
	g++ $(flags) glob.o fpgrowth.o fp_tree.o buffer.o fp_node.o -o fim_all

fim_closed: fpclosed.o fp_tree.o buffer.o fsout.h fp_node.o
	g++ $(flags) glob.o fpclosed.o fp_tree.o buffer.o fp_node.o -o fim_closed

fim_maximal: fpmax.o fp_tree.o buffer.o fsout.h fp_node.o
	g++ $(flags) glob.o fpmax.o fp_tree.o buffer.o fp_node.o -o fim_maximal

fpgrowth.o: fpmax.cpp glob.o fp_tree.h common.h fsout.h
	g++ $(flags) -DFI -c fpmax.cpp -o fpgrowth.o

fpclosed.o: fpmax.cpp glob.o fp_tree.h common.h fsout.h
	g++ $(flags) -DCFI -c fpmax.cpp -o fpclosed.o

fpmax.o: fpmax.cpp glob.o fp_tree.h common.h fsout.h main.h
	g++ $(flags) -DMFI -c fpmax.cpp

glob.o: glob.cpp
	g++ $(flags) -c glob.cpp

fp_tree.o: fp_tree.cpp fp_tree.h common.h fsout.h algos.h
	g++ $(flags) -c fp_tree.cpp

buffer.o: buffer.cpp buffer.h
	g++ $(flags) -c buffer.cpp

fp_node.o: fp_node.cpp fp_node.h fsout.h common.h
	g++ $(flags) -c fp_node.cpp

clean:
	rm -f *.o fim_all fim_closed fim_maximal
