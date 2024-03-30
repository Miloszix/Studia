.EQU wy = 0b11111111
.EQU we = 0b00000000

ldi r16, wy
ldi r17, we

out ddrb, r16
out ddrc, r17
out portc, r16 


led:

sbic pinc, 0
cbi portb, 0
sbis pinc, 0
sbi portb, 0
rjmp led

