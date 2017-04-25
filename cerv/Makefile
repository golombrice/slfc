CC=gcc -O3 
CFLAGS=-I. 
DEPS = 

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libcerv.a: adaptive_model.o arithmetic_decode.o arithmetic_encode.o bit_input.o  bit_output.o  cerv.o 
	ar rsc libcerv.a adaptive_model.o arithmetic_decode.o arithmetic_encode.o bit_input.o  bit_output.o  cerv.o

clean:
	rm *.o *.a
