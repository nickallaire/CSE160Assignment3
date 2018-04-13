default: mmmw 

mm: mm.c
	gcc -O3 mm.c -o mm
msgbench: msgbench.c
	mpicc msgbench.c -o msgbench
mmmw: mmmw.c
	mpicc mmmw.c -o mmmw
clean:
	- /bin/rm mm
	- /bin/rm msgbench
	- /bin/rm mmmw
