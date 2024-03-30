// 2.2

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

	sbis pinb, 0
	sbi porta, 0

	sbic pinb, 0
	cbi porta, 0

stop: rjmp start