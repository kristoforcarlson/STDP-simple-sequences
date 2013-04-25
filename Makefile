all: quads

quads:  Quads_pk_p.c
	gcc -o trips Quads_pk_p.c -lm

clean:
	rm trips Ca* v* glut* I_* p* P*
