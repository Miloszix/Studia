.EQU wy = 0b11111111
.EQU we = 0b00000000

ldi r16, wy
ldi r17, we
ldi r19, 15

out ddrb, r16
out ddrc, r17
out portc, r16 

loop:
ldi r18, 5
	led:
		out portb, r18
		ldi r16, 8
		inner_loop:
			ldi r17,100
			inner_loop2:
				ldi r20,100
				inner_loop3:
					ldi r21,50
					inner_loop4:
					dec r21
					brne inner_loop4
				dec r20
				brne inner_loop3
			dec r17
			brne inner_loop2
		dec r16
		brne inner_loop 
		inc r18
		cp r19, r18
		brne led
rjmp loop