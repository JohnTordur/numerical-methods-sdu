CPPFLAGS      += -Wall -std=c11
PROG = lecture2
$(PROG): Lecture2.cpp
	g++ -o lecture2 Lecture2.cpp -I. $(CFLAGS)
.PHONY : run
run : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG)
