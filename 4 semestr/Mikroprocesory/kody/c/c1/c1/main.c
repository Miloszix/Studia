#define F_CPU 16000000UL
#include <avr/io.h>
#include <avr/delay.h>
void wait_ms(uint16_t milliseconds) {
    while (milliseconds > 0) {
        // Za³aduj wartoœæ pocz¹tkow¹ do Timer1
        TCNT1 = 0xFFF0; // 65520

        // Uruchom Timer1 z preskalerem 1024
        TCCR1B = (1 << 2) | (1 << 0);

        // Poczekaj na ustawienie flagi przepe³nienia
        while (!(TIFR1 & (1 << TOV1)));

        // Zatrzymaj Timer1
        TCCR1B = 0;

        // Wyczyœæ flagê przepe³nienia
        TIFR1 = (1 << TOV1);

        // Zmniejsz licznik milisekund
        milliseconds--;
    }
}
int main()
{
	unsigned int stan = 1;	
	unsigned char val=0;
	DDRB=0xFF; //Set port A as output
	DDRC=0x00;
	PORTC=0xff;
	while(1)
	{
		if ((PINC & (1<<3)) == 0){
			stan = 1;
		}
		if (PINC == 0x7B){
			stan = 2;
		}
		if (PINC == 0x7D){
			stan = 3;
		}
		if (PINC == 0x7E){
			val = 0;
			PORTB = val;
		}
		if (stan==1){
		val++; //Increment value
		PORTB=val; //Send new value to port A
		//_delay_ms(1000);
		wait_ms(1000);
		}
		if (stan==2){
			val--; //Increment value
			PORTB=val; //Send new value to port A
			_delay_ms(1000);
		}
		if (stan==3){
			 //Increment value
			PORTB=val; //Send new value to port A
			//_delay_ms(1000);
		}
	}
	return 0;
}