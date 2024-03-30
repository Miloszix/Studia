// 1.3

begining: 

ldi r20, 0b11111111

out DDRB, r20

sbi portb, 1
sbi portb, 3
sbi portb, 5
sbi portb, 7

stop: rjmp stop
