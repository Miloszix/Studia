.EQU wy = 0b11111111
.EQU we = 0b00000000

ldi r16, wy
ldi r17, we

out ddrb, r16
out ddrc, r17
out portc, r16 

led:

in r18, pinc

eor r18, r16

out portb, r18

rjmp led

