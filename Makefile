CC=g++ -pthread -O3 -std=c++0x
#CFLAGS=-I. -I/usr/local/include/openjpeg-2.2
#CFLAGS=-I. -I/Users/helinp/openjpeg-2.1.2/include/openjpeg-2.1 -L/Users/helinp/openjpeg-2.1.2/lib -Ilars -Llars -I/Users/helinp/lightfield_compression/lossless_journal_march/codegen/lib/lars_wrapper -L/Users/helinp/lightfield_compression/lossless_journal_march/codegen/lib/lars_wrapper
CFLAGS=-I. -I/Users/helinp/openjpeg-2.1.2/include/openjpeg-2.1 -L/Users/helinp/openjpeg-2.1.2/lib 
#CFLAGS=-I. -I/home/helinp/openjpeg/include/openjpeg-2.2 -L/home/helinp/openjpeg/lib 
#CFLAGS=-I. 
DEPS = image.hh avi_reader.hh

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

compressor: main.o adaptive_coder.o bit_output.o bit_input.o golomb_coder.o adaptive_binary_coder.o displacer.o encoder.o coder.o decoder.o
	$(CC) -o compressor main.o adaptive_coder.o bit_output.o bit_input.o golomb_coder.o adaptive_binary_coder.o displacer.o encoder.o coder.o decoder.o cerv/libcerv.a segm/libms.a $(CFLAGS) -lopenjp2 #-llars_wrapper -llars_main 

clean:
	rm *.o
