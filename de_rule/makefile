CPPFLAGS      += -Wall -std=c11
PROG = derule_exc
$(PROG): derule_exc.cpp
	g++ -o derule_exc derule_exc.cpp -I. $(CFLAGS)
.PHONY : run
run : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG)
