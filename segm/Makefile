CC=g++ -O3 
CFLAGS=-I. 
DEPS = 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libms.a: RAList.o  ms.o   msImageProcessor.o msSys.o   msSysPrompt.o  rlist.o 
	ar rsc libms.a RAList.o  ms.o   msImageProcessor.o msSys.o   msSysPrompt.o  rlist.o

clean:
	rm *.o *.a
