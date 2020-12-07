# Template of makefile

# Write below after "TARGET =" the name of your program without ".c".
# (As the name of the sample program is "template.c", the following
# line is "template".)

TARGET = template

# The default compiler is "gcc" with options "-Wall O2".
# You can change the compiler and options by modifying the following
# lines appropriately.

CC= gcc
CFLAGS= -Wall -O2

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).o -lm

$(TARGET).o: $(TARGET).c cpu_time.c
	$(CC) $(CFLAGS) -c $(TARGET).c

clean:
	rm *.o
