// 2.4

ldi		r17, 5
ldi		r19, 0

loop1:
	ldi r18, 6
	loop2:
		dec r18
		CPSE r18,r19
		jmp loop2
	dec r17
	CPSE r17,r19
	jmp loop1