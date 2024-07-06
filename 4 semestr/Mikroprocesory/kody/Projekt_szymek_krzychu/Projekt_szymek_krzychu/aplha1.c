#include <avr/io.h>
#include <util/delay.h>
#include <avr/interrupt.h>

#define F_CPU 16000000UL
#define delay 200

volatile int16_t encoderValue = 0;
volatile uint8_t lastStateCLK;
volatile uint8_t digits[4];
volatile uint8_t position = 0;
volatile uint8_t currentDigit = 0;
volatile uint8_t step = 0;
volatile uint8_t state = 0;
volatile int number[6];
const int code[6] = {1, 2, 3, 4, 5, 6};

void setup(void) {
	DDRC = 0x00;
	PORTC = 0xff;
	DDRB = 0xff;
	DDRD = 0xff;
	PORTD = 0x00;
	
	lastStateCLK = PINC & (1 << 0);  // pocz¹tkowy stan CLK
	
	// Konfiguracja przerwañ
	PCICR |= (1 << PCIE1); // w³¹czenie przerwañ dla pinów C
	PCMSK1 |= (1 << PCINT8) | (1 << PCINT10); // w³¹czenie na pojedynczych pinach
	
	// konfiguracja timera2
	TCCR2A |= (1 << WGM21); // ustawienie w tryb ctc
	TCCR2B |= (1 << CS22) | (1 << CS21); // Ustawienie preskalera na 256
	OCR2A = 61; // 16000000/(256*1000) = 61
	TIMSK2 |= (1 << OCIE2A); // w³¹czenie przerwania przez porównanie
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
	// w³¹czenie odpowiednich segmentów
	PORTD = digit[digit_value];
	// wybranie odpowiedniej pozycji
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
	if (state == 0) {
		// oddzielenie odpowiedniej cyfry
		digits[0] = currentDigit / 1000;
		digits[1] = (currentDigit / 100) % 10;
		digits[2] = (currentDigit / 10) % 10;
		digits[3] = currentDigit % 10;

		show_digit(digits[position], position);

		// zmiana pozycji
		if (position == 3) {
			position = 0;
			} else {
			position++;
		}
		} else if (state == 1) {
		PORTB = 0x00;
		PORTD = digit[5];
		} else if (state == 2) {
		PORTB = 0x00;
		PORTD = digit[1];
	}
}

// odczyt enkodera
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

int main(void) {
	setup();
	sei();
	int lastencoderValue = encoderValue;
	int direction = 1;
	int lastdirection = 1;
	currentDigit = 0;
	for (int i = 0; i < 6; i++) {
		number[i] = 0;
	}

	while (1) {
		while (step < 6) {
			if (lastencoderValue < encoderValue) {
				direction = 1;
			}
			if (lastencoderValue > encoderValue) {
				direction = -1;
			}
			if (lastdirection != direction) {
				step++;
				lastdirection = direction;
				encoderValue = 0;
				lastencoderValue = encoderValue;
				} else {
				if (direction == 1) {
					lastencoderValue = encoderValue;
					number[step] = encoderValue;
					} else {
					lastencoderValue = encoderValue;
					number[step] = -encoderValue;
				}
			}
			if (number[step] >= 100 || number[step] <= -100) {
				encoderValue = 0;
				lastencoderValue = encoderValue;
				number[step] = 0;
			}
			currentDigit = abs(number[step]);
		}

		if (step == 6) {
			uint8_t correct = 1;
			for (int i = 0; i < 6; i++) {
				if (number[i] != code[i]) {
					correct = 0;
					break;
				}
			}
			if (correct) {
				state = 1; // Kod poprawny
				} else {
				state = 2; // Kod niepoprawny
			}
		}
	}
}