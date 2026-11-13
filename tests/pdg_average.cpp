#include "PDGAverage.h"
#include <cmath>
#include <cstring>
#include <new>
#include <stdexcept>

int main() {
    // Poison storage to reproduce the previously uninitialized angular flag.
    alignas(PDGAverage) unsigned char storage[sizeof(PDGAverage)];
    std::memset(storage, 0xff, sizeof(storage));
    auto* average = new (storage) PDGAverage;
    average->setData({dato(.4, .2), dato(1.2, .3)});
    average->CalculateAverage();
    const double expected = (.4/.04 + 1.2/.09) / (1./.04 + 1./.09);
    if (std::abs(average->getAverage()-expected) > 1.e-14)
        throw std::runtime_error("Default PDG averaging must be arithmetic on every rank");
    average->setIsAngle(true);
    average->setData({dato(3.1, .1), dato(-3.1, .1)});
    average->CalculateAverage();
    if (std::abs(std::abs(average->getAverage())-M_PI) > 1.e-14)
        throw std::runtime_error("Explicit circular PDG averaging must remain available");
    average->~PDGAverage();
}
