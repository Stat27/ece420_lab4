CC = mpicc
CFLAGS = -g -Wall

all: main_mpi main_serial datatrim

main_mpi: main_mpi.o Lab4_IO.o
	$(CC) $(CFLAGS) -o main_mpi main_mpi.o Lab4_IO.o -lm

main_serial: main_serial.o Lab4_IO.o
	$(CC) $(CFLAGS) -o main_serial main_serial.o Lab4_IO.o -lm

datatrim: datatrim.o
	gcc -o datatrim datatrim.o

main_mpi.o: main_mpi.c
	$(CC) $(CFLAGS) -c main_mpi.c -o main_mpi.o

main_serial.o: main_serial.c
	$(CC) $(CFLAGS) -c main_serial.c -o main_serial.o

Lab4_IO.o: Lab4_IO.c
	$(CC) $(CFLAGS) -c Lab4_IO.c -o Lab4_IO.o

datatrim.o: datatrim.c
	gcc -c datatrim.c -o datatrim.o

clean:
	rm -f main_mpi main_serial datatrim *.o data_input_* data_output*