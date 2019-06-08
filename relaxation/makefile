CPPFLAGS      += -Wall -std=c11
PROG = relax
$(PROG): relax.cpp
	g++ -o relax relax.cpp -I. $(CFLAGS)
.PHONY : run
run : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG)
.PHONY : plot
plot : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG) | feedgnuplot --domain &
