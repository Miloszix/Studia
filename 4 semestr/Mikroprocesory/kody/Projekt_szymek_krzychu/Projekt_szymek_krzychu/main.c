
#include <avr/io.h>
#include <avr/delay.h>
#include <avr/interrupt.h>

#define F_CPU 16000000UL
#define delay 200

volatile int16_t encoderValue = 0;
volatile uint8_t lastStateCLK;
volatile uint8_t digits[4];
volatile uint8_t position=0;
volatile uint8_t number;

void setup(void) {
	DDRC = 0x00;
    PORTC = 0xff;
	DDRB = 0xff;
	DDRD = 0xff;
	PORTD = 0x00;
    
    lastStateCLK = PINC & (1 << 0);  //pocz¹tkowy stan CLK
    
    // Konfiguracja przerwañ
    PCICR |= (1 << PCIE1);
    PCMSK1 |= (1 << PCINT8) | (1 << PCINT9) | (1 << PCINT10);
	
    TCCR2A |= (1 << WGM21);
    TCCR2B |= (1 << CS22) | (1 << CS21); // Ustawienie preskalera na 256
    OCR2A = 61; //16000000/(256*1000) = 61
    TIMSK2 |= (1 << OCIE2A);
}

const uint8_t digit[10] = {
	0b11000000, // 0
    0b11111001, // 1
    0b10100100, // 2
    0b10110000, // 3
    0b10011001, // 4
    0b10010010, // 5
    0b10000010, // 6
    0b11111000, // 7
    0b10000000, // 8
    0b10010000  // 9
};

void show_digit(uint8_t digit_value, uint8_t position) {
	// Wy³¹czenie cyfr
	PORTB = 0xff;
	// ustawienie segmentów wed³ug powy¿szej konfiguracji
	PORTD = digit[digit_value];
	// selekcja odpowiedniej cyfry
	switch (position) {
		case 0:
		PORTB = ~(1 << 0);
		break;
		
		case 1:
		PORTB = ~(1 << 1);
		break;
		
		case 2:
		PORTB = ~(1 << 2);
		break;
		
		case 3:
		PORTB &= ~(1 << 3);
		break;
	}
}

ISR(TIMER2_COMPA_vect) {
		//reszty z dzielenia wyœwietlane jako cyfry
		digits[0] = number / 1000;
		digits[1] = (number / 100) % 10;
		digits[2] = (number / 10) % 10;
		digits[3] = number % 10;

		show_digit(digits[position], position);
		
	if (position==3){
		position=0;
		}else{
		position++;
	}
}

ISR(PCINT1_vect) {
	// Odczytaj bie¿¹cy stan CLK
	uint8_t currentStateCLK = PINC & (1 << 0);
	_delay_us(delay);
	
		if (currentStateCLK != lastStateCLK) {  
			uint8_t stateDT = PINC & (1 << 1); 
			if (currentStateCLK == 0) {
				if (stateDT == 0) {
					encoderValue++; // Obrót ze wskazówkami zegara
				} else {
					encoderValue--; // Obrót w przeciwnym kierunku do wskazówek zegara
				}
			}
		}
	lastStateCLK = currentStateCLK;
}


int main(void)
{
	setup();
	sei();
		
    while (1) 
    {
		number=encoderValue;
    }
}

