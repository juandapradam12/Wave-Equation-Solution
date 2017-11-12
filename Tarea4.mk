# Tarea4.mk 

# compilo con make -f Tarea4.mk

Resultados_hw4.pdf : Resultados_hw4.tex *.png
	pdflatex $<

*.png : *.dat Plots.py
	python Plots.py

*.dat : Ondas.x 
	./Ondas.x 

Ondas.x : Ondas.c
	cc Ondas.c -o Ondas.x -lm

clean :
	rm *.dat  *.x *.png *.log *.aux 
