// 1.1

.EQU decimal=10
.EQU hex=0x99
.EQU binary=0b01001
.EQU ascii='A'; 

ldi r16, decimal
ldi r17, hex
ldi r18, binary 
ldi r19, ascii

stop: rjmp stop
	
