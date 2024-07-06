/*
 * Projekt.c
 *
 * Created: 27.05.2024 18:50:56
 * Author : Mi³osz Mynarczuk, Dawid Makowski
 */
 
#define F_CPU 16000000UL
#include <avr/io.h>
#include <avr/interrupt.h>

#define DEBOUNCE_TIME 1 //ms

volatile int16_t target; // mno¿nik encoder value
volatile int16_t steps; // po³o¿enie silnika
volatile int16_t encoderValue = 0;
volatile uint8_t lastStateCLK;
volatile uint8_t state = 0;
volatile uint8_t position=0;

void setup(void) {
	/* 
	PC0 clk PC1 DT PC2 SW PC3 D1 PC4 D2 PC5 D3
	PE0 S1 PE1 S2 PE2 S3 PE3 S4
	D0-6 segmenty PD7 D4
	PB silnik
    Ustawienie portów C0-2 i E jako wejœæ, B i D C3-8 jako wyjœæ 
    */
    
	DDRC = 0b11111000;
    PORTC = 0b00000111;	
	DDRB = 0xff;
	DDRD = 0xff;
	DDRE = 0x00;
	PORTE = 0xff;
	PORTD = 0x00;
    
    lastStateCLK = PINC & (1 << 0);  //pocz¹tkowy stan CLK
    
    // Konfiguracja przerwañ
    PCICR |= (1 << PCIE1) | (1 << PCIE3); // W³¹czenie przerwañ dla pinów C i E
    PCMSK1 |= (1 << PCINT8) | (1 << PCINT9) | (1 << PCINT10); // W³¹czenie przerwañ dla PC0-2
	PCMSK3 |= (1 << PCINT24) | (1 << PCINT25) | (1 << PCINT26) | (1 << PCINT27); // W³¹czenie przerwañ dla PE0-3
	
    // Ustawienie Timer0 w tryb CTC
    TCCR0A |= (1 << WGM01);
    TCCR0B |= (1 << CS01) | (1 << CS00); // Ustawienie preskalera na 256
    OCR0A = 249; //16000000/(64*1000) = 249 ustawienie wartoœci rejestru porównania
    TIMSK0 |= (1 << OCIE0A); //w³¹cza przerwanie porównania
}

const uint8_t digit[11] = {
	  ~((1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5)),          // 0
	  ~((1 << 1) | (1 << 2)),                                                      // 1
	  ~((1 << 0) | (1 << 1) | (1 << 3) | (1 << 4) | (1 << 6)),                     // 2
	  ~((1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 6)),                     // 3
	  ~((1 << 1) | (1 << 2) | (1 << 5) | (1 << 6)),                                // 4
	  ~((1 << 0) | (1 << 2) | (1 << 3) | (1 << 5) | (1 << 6)),                     // 5
	  ~((1 << 0) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6)),          // 6
	  ~((1 << 0) | (1 << 1) | (1 << 2)),                                           // 7
	  ~((1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4) | (1 << 5) | (1 << 6)), // 8
	  ~((1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 5) | (1 << 6)), //9
	  ~(1 << 6) // -
};


void wait_ms(uint16_t milliseconds) {
	while (milliseconds > 0) {
		// Ustaw Timer1 w tryb CTC
		TCCR1A = 0; // Wyzerowanie rejestru kontrolnego A
		TCCR1B = (1 << WGM12) | (1 << CS11) | (1 << CS10); // WGM12: CTC tryb, CS11 | CS10: Preskaler 64

		// Ustaw wartoœæ porównania dla OCR1A, aby uzyskaæ 1 ms
		OCR1A = 249; // Liczba cykli na 1 ms - 1

		// Wyczyœæ licznik
		TCNT1 = 0;

		// Poczekaj na ustawienie flagi porównania
		while (!(TIFR1 & (1 << OCF1A)));

		// Wyczyœæ flagê porównania
		TIFR1 = (1 << OCF1A);

		// Zmniejsz licznik milisekund
		milliseconds--;
	}
}

void display_digit(uint8_t digit_value, uint8_t position) {
	// Wy³¹czenie cyfr
	PORTD |= (1 << 7);
	PORTC = ((1 << 3) | (1 << 4) | (1 << 5));
	// ustawienie segmentów wed³ug powy¿szej konfiguracji
	PORTD = digit[digit_value];
	// selekcja odpowiedniej cyfry
	switch (position) {
		case 0:
		PORTC = ~(1 << 3);
		break;
		
		case 1:
		PORTC = ~(1 << 4);
		break;
		
		case 2:
		PORTC = ~(1 << 5);
		break;
		
		case 3:
		PORTD &= ~(1 << 7);
		break;
	}
}

void display_number(int16_t number) {
	uint8_t digits[4];
	//warunek dla liczb dodatnich
	if (number>=0){                   
		digits[0] = number / 1000;	
	}
	//wyœwietlenie minusa dla ujemnych
	if (number<0){      
		digits[0] = 10;             
		number = -number;
	}
	//reszty z dzielenia wyœwietlane jako cyfry
	digits[1] = (number / 100) % 10; 
	digits[2] = (number / 10) % 10;
	digits[3] = number % 10;

	display_digit(digits[position], position);
	
}

//przerwania na poruszenie enkodera
ISR(PCINT1_vect) {
    // Odczytaj bie¿¹cy stan CLK
    uint8_t currentStateCLK = PINC & (1 << 0);
	wait_ms(DEBOUNCE_TIME); // zapobieganie glitchom podczas przekrêcania
	
	//przypadek dla programu 0 kiedy switch na enkoderze jest wciœniêty	
	if (((PINC&(1<<2))==0)&(state==0)){
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni³o siê CLK
			uint8_t stateDT = PINC & (1 << 1); // Pobranie wartoœci DT
			// Okreœlenie kierunku obrotu enkodera CLK zmienia siê z 1 na 0
			if (currentStateCLK == 0) {  
				if (stateDT == 0) {
					encoderValue=encoderValue+10; // Obrót ze wskazówkami zegara
				} else {
					encoderValue=encoderValue-10; // Obrót w przeciwnym kierunku do wskazówek zegara
				}
			}
			target=encoderValue*5;
		}
	}else{  
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni³o siê CLK
			uint8_t stateDT = PINC & (1 << 1); // Pobranie wartoœci DT
			// Okreœlenie kierunku obrotu enkodera CLK zmienia siê z 1 na 0     
			if (currentStateCLK == 0) {  
				if (stateDT == 0) {
					encoderValue++; // Obrót ze wskazówkami zegara
				} else {
					encoderValue--; // Obrót w przeciwnym kierunku do wskazówek zegara
				}
			}
			target=encoderValue*5;
		}
	}
	// ograniczenie wartoœci enkoder dla programu 1 i 2, od 0 do 9
	if((state==1)|(state==2)){
		if (encoderValue>9){
				encoderValue=9;
			}
			if (encoderValue<0){
				encoderValue=0;
			}
	}
    lastStateCLK = currentStateCLK; // Aktualizacja ostatniego stanu CLK
}

//przerwania na wciœniecie przycisku
ISR(PCINT3_vect){
	if ((PINE & (1 << 0))==0){
		state = 0;
		encoderValue =0;
		steps = 0;
		target = 0;
	}
	if ((PINE & (1 << 1))==0){
		state = 1;
		encoderValue =0;
	}
	if ((PINE & (1 << 2))==0){
		state = 2;
		encoderValue =0;
	}
	if ((PINE & (1 << 3))==0){
		state = 3;
	}
}

//przerwanie na wyœwietlenie pojedynczej cyfry
ISR(TIMER0_COMPA_vect) {
	display_number(encoderValue);
	if (position==3){
		position=0;
	}else{
		position++;
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
		
		case 5: // od³¹czenie silnika
		PORTB =0x00;
		break;
	}
}

int main(void) {

    setup(); // Ustawienie portów i przerwañ
	    
    sei(); // W³¹cz przerwania globalne
	
	int i= 0;
	int speed;
    
    while (1) {
		// obrót do stanu enkodera z minimalnym opóŸnieniem
		if(state==0){		
			while(target>steps){
				steps++;
				if(i==4){
					i=1;
				}else{
					i++;
				}
				stepMotor(i);
				wait_ms(2);
			}
			while(target<steps){
				steps--;
				if(i<=1){
					i=4;
				}else{
					i--;
				}
				stepMotor(i);
				wait_ms(2);	
			}
		}
		
		//sta³y obrót ze wskazówkami zegara
		if(state==1){			
			speed = 10 - encoderValue;
			stepMotor(i);
			
			if(i==4){
				i=1;
			}else{
				i++;
			}		
			//opó¿nienie do kontroli prêdkoœci
			for(uint8_t t = 0; t<speed; t++ ){
				wait_ms(2);
			}
		}
		
		//sta³³y obrót w kierunku przeciwnym do wskazówek zegara
		if(state==2){	
			speed = 10 - encoderValue;
			stepMotor(i);
			
			if(i<=1){
				i=4;
			}else{
				i--;
			}
			//opó¿nienie do kontroli prêdkoœci
			for(uint8_t t = 0; t<speed; t++ ){
				wait_ms(2);
			}
		}
		
		//zatrzymanie programów i od³¹czenie silnika 
		if(state==3){
			encoderValue=0;
			stepMotor(5);
		}		
	}	
    return 0;
}