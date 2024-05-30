#define F_CPU 16000000UL
#include <avr/io.h>
#include <avr/delay.h>

int main()
{
	unsigned int stan = 1;	
	unsigned char val=0;
	DDRB=0xFF; //Set port A as output
	DDRC=0x00;
	PORTC=0xff;
	while(1)
	{
		if (PINC == 0x77){
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
		while (stan==1){
		val++; //Increment value
		PORTB=val; //Send new value to port A
		_delay_ms(1000);
		}
		while (stan==2){
			val--; //Increment value
			PORTB=val; //Send new value to port A
			_delay_ms(1000);
		}
		while (stan==3){
			 //Increment value
			PORTB=val; //Send new value to port A
			_delay_ms(1000);
		}
	}
	return 0;
}