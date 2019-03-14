CPPFLAGS      += -Wall -std=c11
PROG = lecture3
$(PROG): Lecture3.cpp
	g++ -o lecture3 Lecture3.cpp -I. $(CFLAGS)
.PHONY : run
run : $(PROG) # assuming $(PROG) is the name of your program
	./$(PROG)
