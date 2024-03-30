#include <iostream>

class Point {
private:
    int x;
    int y;
public:
    Point() {}

    Point(int x1, int y1) {
        x = x1;
        y = y1;
    }

    Point(const Point & other) {
        x = other.x;
        y = other.y;
    }

    void setX(int x1) {
        x = x1;
    }

    void setY(int y1) {
        y = y1;
    }

    int getX() const {
        return x;
    }

    int getY() const {
        return y;
    }

    bool operator==(const Point & other) const {
        return x == other.x && y == other.y;
    }

    void move(int a, int b) {
        x += a;
        y += b;
    }
};

int main() {
    Point Milosz;

    Point Mynarczuk(5, 10);

    Point maj16 = Mynarczuk;

    Milosz.setX(3);

    Milosz.setY(7);

    std::cout << "Wspolrzedne punktu p1: (" << Milosz.getX() << ", " << Milosz.getY() << ")" << std::endl;
    std::cout << "Wspolrzedne punktu p2: (" << Mynarczuk.getX() << ", " << Mynarczuk.getY() << ")" << std::endl;
    std::cout << "Wspolrzedne punktu p3: (" << maj16.getX() << ", " <<maj16.getY() << ")" << std::endl;


    if(Milosz == Mynarczuk) {
        std::cout << "Punkt p1 jest równy punktowi p2" << std::endl;
    } else {
        std::cout << "Punkt p1 nie jest równy punktowi p2" << std::endl;
    }

    
    Mynarczuk.move(2, -5);
    Milosz.move(1, 5);
    maj16.move(-3, 2);

    std::cout << "Wspolrzedne punktu p1: (" << Milosz.getX() << ", " << Milosz.getY() << ")" << std::endl;
    std::cout << "Wspolrzedne punktu p2: (" << Mynarczuk.getX() << ", " << Mynarczuk.getY() << ")" << std::endl;
    std::cout << "Wspolrzedne punktu p3: (" << maj16.getX() << ", " <<maj16.getY() << ")" << std::endl;

}
