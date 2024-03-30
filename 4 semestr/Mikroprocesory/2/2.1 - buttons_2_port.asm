// 2.1

ldi r16, $ff;
out PORTa, r16;
ldi r16, $ff;
out DDRa, r16;
ldi r16, $00;
out PORTa, r16;


ldi r16, $ff;
out PORTb, r16;
ldi r16, $00;
out DDRb, r16;
in r16,PINb;

start:

	in r16, PINb;
	COM r16;
	out PORTa, r16; 

stop: rjmp start