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

#define DEBOUNCE_TIME 100 //us

volatile int16_t target; // mno�nik encoder value
volatile int16_t steps; // po�o�enie silnika
volatile int16_t encoderValue = 0;
volatile uint8_t lastStateCLK;
volatile uint8_t state = 0;

void setup(void) {
	/* 
	PC0 clk PC1 DT PC2 SW PC3 D1 PC4 D2 PC5 D3
	PE0 S1 PE1 S2 PE2 S3 PE3 S4
	D0-6 segmenty PD7 D4
	PB silnik
    Ustawienie port�w C0-2 i E jako wej��, B i D C3-8 jako wyj�� 
    */
    
	DDRC = 0b11111000;
    PORTC = 0b00000111;	
	DDRB = 0xff;
	DDRD = 0xff;
	DDRE = 0x00;
	PORTE = 0xff;
	PORTD = 0x00;
    
    lastStateCLK = PINC & (1 << 0);  //pocz�tkowy stan CLK
    
    // Konfiguracja przerwa�
    PCICR |= (1 << PCIE1) | (1 << PCIE3); // W��czenie przerwa� dla pin�w C i E
    PCMSK1 |= (1 << PCINT8) | (1 << PCINT9) | (1 << PCINT10); // W��czenie przerwa� dla PC0-2
	PCMSK3 |= (1 << PCINT24) | (1 << PCINT25) | (1 << PCINT26) | (1 << PCINT27); // W��czenie przerwa� dla PE0-3
	
	
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


void display_digit(uint8_t digit_value, uint8_t position) {
	// Wy��cz wszystkie cyfry
	PORTD |= (1 << 7);
	PORTC = ((1 << 3) | (1 << 4) | (1 << 5));

	// Ustaw segmenty na odpowiedni� cyfr�
	PORTD = digit[digit_value];

	// W��cz odpowiedni� cyfr�
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
	//wy�wietlenie minusa dla ujemnych
	if (number<0){      
		digits[0] = 10;             
		number = -number;
	}
	//reszty z dzielenia wy�wietlane jako cyfry
	digits[1] = (number / 100) % 10; 
	digits[2] = (number / 10) % 10;
	digits[3] = number % 10;

		for (uint8_t n = 0; n < 4; n++) {
			display_digit(digits[n], n);
			_delay_ms(1); // Kr�tki czas op�nienia dla multipleksowania
		}
}


ISR(PCINT1_vect) {
    // Odczytaj bie��cy stan CLK
    uint8_t currentStateCLK = PINC & (1 << 0);
	
	_delay_us(DEBOUNCE_TIME); // zapobieganie glitchom podczas przekr�cania
	
	if (((PINC&(1<<2))==0)&state==0){ //Sprawdzenie czy switch jest wci�ni�ty i jest to program zerowy i wtedy zmiana jest o 10 a nie 1
		
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni�o si� CLK
		
			uint8_t stateDT = PINC & (1 << 1); // Pobranie warto�ci DT
        
			if (currentStateCLK == 0) {  // Okre�lenie kierunku obrotu enkodera CLK zmienia si� z 1 na 0
				if (stateDT == 0) {
					encoderValue=encoderValue+10; // Obr�t ze wskaz�wkami zegara
				} else {
					encoderValue=encoderValue-10; // Obr�t w przeciwnym kierunku do wskaz�wek zegara
				}
			}
			target=encoderValue;
		}
	}else{
    
		if (currentStateCLK != lastStateCLK) {     // Sprawdzenie czy zmieni�o si� CLK
		
			uint8_t stateDT = PINC & (1 << 1); // Pobranie warto�ci DT
        
			if (currentStateCLK == 0) {  // Okre�lenie kierunku obrotu enkodera CLK zmienia si� z 1 na 0
				if (stateDT == 0) {
					encoderValue++; // Obr�t ze wskaz�wkami zegara
				} else {
					encoderValue--; // Obr�t w przeciwnym kierunku do wskaz�wek zegara
				}
			}
			target=encoderValue;
		}
	}
	// ograniczenie warto�ci enkoder dla programu 1 i 2, od 0 do 9
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
		case 5: // od��czenie silnika
		PORTB =0x00;
		break;
	}
}


int main(void) {

    setup(); // Inicjalizacja enkodera
	    
    sei(); // W��cz przerwania globalne
	int i= 0;
	int speed;
    
    while (1) {
		// obr�t do stanu enkodera
		if(state==0){		
			display_number(encoderValue);
			while(target>steps){
				steps++;
				if(i==4){
					i=1;
				}else{
					i++;
				}
				stepMotor(i);
				_delay_ms(2);
			}
			while(target<steps){
				steps--;
				if(i<=1){
					i=4;
				}else{
					i--;
				}
				stepMotor(i);
				_delay_ms(2);	
			}
		}
		//sta�y obr�t ze wskaz�wkami zegara
		if(state==1){

			display_digit(encoderValue,3); //wy�wietlenie pr�dko�ci
			speed = 10 - encoderValue;
			stepMotor(i);
			
			if(i==4){
				i=1;
			}else{
				i++;
			}
			//op�nienie do kontroli pr�dko�ci
			for(uint8_t t = 0; t<speed; t++ ){
				_delay_ms(2);
			}
		}
		//sta��y obr�t w kierunku przeciwnym do wskaz�wek zegara
		if(state==2){	
			
			display_digit(encoderValue,3); //wy�wietlenie pr�dko�ci
			speed = 10 - encoderValue;
			stepMotor(i);
			
			if(i<=1){
				i=4;
			}else{
				i--;
			}
			//op�nienie do kontroli pr�dko�ci
			for(uint8_t t = 0; t<speed; t++ ){
				_delay_ms(2);
			}
		}
		//zatrzymanie program�w i od��czenie silnika 
		if(state==3){
			display_number(0000);
			stepMotor(5);
		}		
	}
		
	
    return 0;
}