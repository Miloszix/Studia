.EQU ddr=0b11111111
.EQU port=0b01010101

ldi r16, ddr
ldi r17, port

out ddrb, r16
out portb, r17

stop: rjmp stop