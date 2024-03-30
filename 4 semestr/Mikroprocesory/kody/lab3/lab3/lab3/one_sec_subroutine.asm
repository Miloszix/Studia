ldi r16, high(ramend) ;load sph
out sph, r16
ldi r16, low(ramend) ;load spl
out spl, r16

ldi r17, $ff
out ddrb, r17


start:
call wait
inc r16
out portb, r16
rjmp start



wait:
	push r16
	ldi r16,250
		loop1:
			push r17
			ldi r17, 250
				loop2:
					push r18
					ldi r18, 40
						outter_loop:
							dec r18
						brne outter_loop
					pop r18
					dec r17
				brne loop2
			pop r17
			dec r16
		brne loop1
	pop r16
ret