#include <avr/io.h>
/*
*Implement wait function
*/
void wait(unsigned int delay)
{
	volatile unsigned int i; //volatile prevents from i optimization
	for(i=0; i<delay; i++); //empty loop as delay
}
int main()
{
	unsigned char val=0;
	DDRB=0xFF; //Set port A as output
	while(1)
	{
		val++; //Increment value
		PORTB=val; //Send new value to port A
		wait(5000);
	}
	return 0;
}