CPPFLAGS      += -Wall -std=c11
PROG = shooting
$(PROG): shooting.cpp
	g++ -o shooting shooting.cpp -I. $(CFLAGS)
.PHONY : run
run : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG)
.PHONY : plot
plot : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG) | feedgnuplot --domain &
