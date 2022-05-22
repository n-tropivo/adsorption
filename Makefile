all:
	g++ IAS.cpp -o IAS.x -lgsl

clean:
	rm *.x