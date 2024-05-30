/*
 * Projekt.c
 *
 * Created: 27.05.2024 18:50:56
 * Author : milos
 */
 
#define F_CPU 16000000UL
#include <avr/io.h>
#include <avr/interrupt.h>
#include <avr/delay.h>

#define DEBOUNCE_TIME 1

volatile int16_t target; // po³o¿enie do którego d
volatile int16_t steps; // po³o¿enie silnika
volatile int16_t encoderValue = 0;
volatile uint8_t lastStateCLK;
int state = 0;

void setup(void) {
	/*C enkoder
	E pzyciski
	D ledy
	B silnik
    Ustawienie portów C i E jako wejœæ, B i D jako wyjœæ */
    
	DDRC = 0x00;
    PORTC = 0xff;	
	DDRB = 0xff;
	DDRD = 0xff;
	DDRE = 0x00;
	PORTE = 0xff;
    
    lastStateCLK = PINC & (1 << 0);  //pocz¹tkowy stan CLK
    
    // Konfiguracja przerwañ
    PCICR |= (1 << PCIE1) | (1 << PCIE3); // W³¹czenie przerwañ dla pinów C i E
    PCMSK1 |= (1 << PCINT8) | (1 << PCINT9) | (1 << PCINT10); // W³¹czenie przerwañ dla PC0-2
	PCMSK3 |= (1 << PCINT24) | (1 << PCINT25) | (1 << PCINT26) | (1 << PCINT27); // W³¹czenie przerwañ dla PE0-3
}

ISR(PCINT1_vect) {
    // Odczytaj bie¿¹cy stan CLK
    uint8_t currentStateCLK = PINC & (1 << 0);
	
	_delay_ms(DEBOUNCE_TIME); // zapobieganie glitchom podczas przekrêcania
	
	if ((PINC&(1<<2))==0){ //Sprawdzenie czy switch jest wciœniêty i wtedy zmiana jest o 10 a nie 1
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni³o siê CLK
		
			uint8_t stateDT = PINC & (1 << 1); // Pobranie wartoœci DT
        
			if (currentStateCLK == 0) {  // Okreœl kierunek obrotu enkodera CLK zmienia siê z 1 na 0
				if (stateDT == 0) {
					encoderValue=encoderValue+10; // Obrót ze wskazówkami zegara
				} else {
					encoderValue=encoderValue-10; // Obrót w przeciwnym kierunku do wskazówek zegara
				}
			}
			PORTD=encoderValue;	
			target=encoderValue*5;
		}
	}else{
    
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni³o siê CLK
		
			uint8_t stateDT = PINC & (1 << 1); // Pobranie wartoœci DT
        
			if (currentStateCLK == 0) {  // Okreœl kierunek obrotu enkodera CLK zmienia siê z 1 na 0
				if (stateDT == 0) {
					encoderValue++; // Obrót ze wskazówkami zegara
				} else {
					encoderValue--; // Obrót w przeciwnym kierunku do wskazówek zegara
				}
			}
			PORTD=encoderValue;
			target=encoderValue*5;
		}
	}    
    lastStateCLK = currentStateCLK; // Aktualizacja ostatniego stanu CLK
}

ISR(PCINT3_vect){
	if ((PINE & (1 << 0))==0){
		state = 0;
		encoderValue =0;
		steps = 0;
		target = 0;
	}
	if ((PINE & (1 << 1))==0){
		state = 1;
	}
	if ((PINE & (1 << 2))==0){
		state = 2;
	}
	if ((PINE & (1 << 3))==0){
		state = 3;
	}
}


void stepMotor(int step) {
	switch (step) {
		case 1:
		PORTB = (1 << 0);
		break;
		case 2:
		PORTB = (1 << 1);
		break;
		case 3:
		PORTB = (1 << 2);
		break;
		case 4:
		PORTB = (1 << 3);
		break;
	}
}


int main(void) {

    setup(); // Inicjalizacja enkodera
	    
    sei(); // W³¹cz przerwania globalne
	int i= 1;
    
    while (1) {
		if(state==0){		
			while(target>steps){
				stepMotor(i);
				steps++;
				if(i==4){
					i=1;
				}else{
					i++;
				}
				_delay_ms(2);
			}
			while(target<steps){
				stepMotor(i);
				steps--;
				if(i==1){
					i=4;
				}else{
					i--;
				}
				_delay_ms(2);	
			}
		}	
	}
    return 0;
}