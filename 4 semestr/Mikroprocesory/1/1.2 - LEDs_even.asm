// 1.2

 begining: 

ldi r20, 0b11111111
ldi r21, 0b10101010

out DDRB, r20
out PORTB, r21

stop: rjmp stop
