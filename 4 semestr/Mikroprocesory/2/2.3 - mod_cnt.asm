// 2.3

ldi		r17,$ff
out		ddra,r17
ldi		r20, 15
ldi		r18, 5

loop1: 
	out porta, r18
	inc r18

	ldi r19, 100
	loop2: 
	    ldi r16,50
		outter_loop:
		    ldi r17,50
			inner_loop: nop
				dec r17
			brne inner_loop
		dec r16 
		brne outter_loop
	dec r19
	brne loop2
	cp r18, r20
brne loop1