.DSEG ;start data memory
.ORG 0x100 ;beginning of the SRAM memory in data space
var1: .BYTE 1 ; reserve 1 byte to var1
.CSEG ;program memory again

ldi r16, high(ramend) ;load sph
out sph, r16
ldi r16, low(ramend) ;load spl
out spl, r16

ldi r17, $ff
out ddrb, r17

ldi r16, $00
ldi r18, 5
sts var1, r18

prg:
call wait
inc r16
out portb, r16
rjmp prg



wait:
	push r16
	push r17
	lds r19, var1
	loop0:
	ldi r16,200
		loop1:
			ldi r17, 32
				loop2:
					ldi r18, 250
						outter_loop:
							nop
							nop
							dec r18
						brne outter_loop
					dec r17
				brne loop2
			dec r16
		brne loop1
	dec r19
	brne loop0
	pop r17
	pop r16
ret