// 5.3

.org 0

ldi R16, HIGH(RAMEND)
out SPH, R16
ldi R16, LOW(RAMEND)
out SPL, R16

;inicjalizacja portu a i c jako wyjsc
ldi r17, 0xFF
out ddrb, r17
out ddrd, r17
out ddrc, r17
jmp start1
;inicjalizacja tablicy przechowujacej kolejno liczby 0 1 2... F na wyswietlacz 7-segment
.org 0x32
prime: .DB ~0x7e, ~0x30, ~0x6d, ~0x79, ~0x33, ~0x5b, ~0x5f, ~0x70, ~0x7F, ~0x7B, ~0x77, ~0x1f, ~0x4e, ~0x3d, ~0x4f, ~0x47, ~0x47
.org 0x80

start1:
;inicjalizacja wskaznika (rejestru z) na tablice
ldi zl, low(2*prime);
ldi zh, high(2*prime);
ldi r20, 0x38 ; mniej znaczacy bit temp w stopniach fahrenheita w hexie
ldi r21, 0xff ; bardziej znaczacy bit temp w stopniach fahrenheita w hexie
subi r20,0x20 ; Subtract low byte
sbci r21,0x00 ; Subtract with carry high byte
brvs over1
over2:
brmi sign1 ;jesli wynik odejmowania ujemny to przejdz do sign1, jesli dodatni to rzejdz dalej
ldi r17,0x00
out portb, r17
sign2:
;ldi r17, 0x7f;
;and r21, r17;
mov r30, r20
mov r31, r21
mov r29, r21

ldi r17, 0x01
and r29, r17
swap r29
lsl r29
lsl r29
lsl r29
ldi r17, 0x7f
or r29, r17
lsr r20
ldi r17, 0x80
or r20, r17
and r20, r29
asr r21
mov r29, r31
ldi r17, 0x0f
and r29, r17
swap r29

lsr r30
lsr r30
lsr r30
lsr r30
asr r31
asr r31
asr r31
asr r31
ldi r17, 0xf0
or r30, r17
ldi r17, 0x0f
or r29, r17
and r30, r29
add r20, r30
adc r21, r31

l1: ;nasz licznik
call display
rjmp l1

over1:
ldi r17, 0xFF
out portd, r17
jmp over2
sign1:
ldi r17, 0x01
out portb, r17 ; zaswiec diode na pc0 zeby pokazac ujemny znak
neg r20
neg r21
subi r21, 0x01
jmp sign2

wait:	;ustawione na 1/20s
ldi R16, 2
loop1:
ldi R17, 50
loop2: 
ldi R18, 100
loop3:
dec R18
brne loop3
nop
dec R17
brne loop2
nop
dec R16
brne loop1
nop
ret

display:
ldi r26,5 ; r26 sluzy aby cyfry zminialy sie co 1/4s - zwykla petla
loop0:
call seg1 
call wait
call seg2
call wait
call seg3
call wait
call seg4
call wait
dec r26
brne loop0
ret

seg1:
ldi zl, low(2*prime) 
ldi zh, high(2*prime);ustawiamy wskaznik na zerowy element tablicy, pierwszy znak znajduje sie pod pierwszym elementem talbicy
ldi r25,~0x01 ;wybieramy wyswietlacz pierwszy pod portem pc0
out portc, r25
mov r23, r20 ;kopiujemy rejestr r20 (licznik) aby nei zgubic jego wartosci - mozna zrzucic na stos rownie dobrze
andi r23, 0x0f ;wybieramy tylko 4 najnizsze bity 
ldi r27, 0x01 ;te dwie komendy wlasnie maja nam pomoc przesunac sie o jeden element dalej w tablicy tzn. jesli r20 bylo 0 to przesuwamy sie
add r23, r27 ;na pierwszy element tablicy czyli 0x7e, jesli r20 bylo 1 to dodajemy 1 i przesuwamy sie na drugi element tablicy czyli 0x30
l2:
lpm r24, z+ 
dec r23
brne l2
out portd, r24 ;wyswietlamy dany element tablicy
ret

seg2:
ldi zl, low(2*prime);
ldi zh, high(2*prime);
ldi r25,~0x02
out portc, r25
mov r23, r20
andi r23, 0xf0
swap r23
ldi r27, 0x01
add r23, r27
l3:
lpm r24, z+
dec r23
brne l3
out portd, r24
ret
seg3:
ldi zl, low(2*prime);
ldi zh, high(2*prime);
ldi r25,~0x04
out portc, r25
mov r23, r21
andi r23, 0x0f
ldi r27, 0x01
add r23, r27
l5:
lpm r24, z+
dec r23
brne l5
out portd, r24
ret
seg4:
ldi zl, low(2*prime);
ldi zh, high(2*prime);
ldi r25,~0x08
out portc, r25	
mov r23, r21
andi r23, 0xf0
swap r23
ldi r27, 0x01
add r23, r27
l7:
lpm r24, z+
dec r23
brne l7
out portd, r24
ret


start:
jmp start