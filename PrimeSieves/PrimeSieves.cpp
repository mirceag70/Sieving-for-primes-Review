#include "Helper.h"

uint64_t Trial_Division(const uint64_t Nmax, uint64_t vPrimes[], void*, void*);
uint64_t Eratosthenes_Basic(const uint64_t Nmax, bool vPrimes[], void*, void*);
uint64_t Singleton_optimized311(uint64_t Nmax, uint64_t IP[], uint64_t IQ[], int32_t JQ[]);
uint64_t Mairson_Standard(uint64_t Nmax, uint64_t rlink[], uint64_t llink[], uint64_t DELETE[]);
uint64_t GalePratt(uint64_t Nmax, uint8_t C[], void*, void*);
uint64_t GriesMisra_Standard(uint64_t Nmax, uint64_t rlink[], uint64_t llink[], void*);
uint64_t SoP(const uint64_t Nmax, uint64_t vPrimes[], void*, void*);
uint64_t Atkin(const uint64_t Nmax, bool vPrimes[], void*, void*);


void Try357(const uint64_t Nmax);
void TryBen(const uint64_t Nmax);

constexpr auto LIMIT = 1'000'000'000;

int main()
{
    std::locale mylocale("");   // get global locale 
    std::cout.imbue(mylocale);  // imbue global locale for thousands delimiter


    //Try_Sieve<uint64_t, int, int, LIMIT / 4>    // for LIMIT >= 100 !!!
    //    (LIMIT, "Trial Division", &Trial_Division, true);
    //Try_Sieve<bool, int, int, LIMIT + 1>
    //    (LIMIT, "Basic Eratosthenes", &Eratosthenes_Basic, true);
    //Try_Sieve<uint64_t, uint64_t, int32_t, LIMIT / 4, LIMIT / 4, LIMIT / 4>
    //    (LIMIT, "Chartres G2 - Singleton 356", &Singleton_optimized311, true);
    //Try_Sieve<uint64_t, uint64_t, uint64_t, LIMIT + 1, LIMIT + 1, LIMIT + 1>
    //    (LIMIT, "Mairson Standard", & Mairson_Standard, true);
    //Try_Sieve<uint8_t, int, int, LIMIT + 1>
    //    (LIMIT, "Gale & Pratt", &GalePratt, true);
    //Try_Sieve<uint64_t, uint64_t, int, LIMIT + 1, LIMIT + 1>
    //    (LIMIT, "Gries/Misra Standard", &GriesMisra_Standard, true);
    //Try_Sieve<uint64_t, int, int, LIMIT + 1000>
    //    (LIMIT, "Pritchard", &SoP, true);
    //Try_Sieve<bool, int, int, LIMIT>
    //    (LIMIT, "Basic Atkin", &Atkin, true);
    //TryBen(LIMIT);
    Try357(LIMIT);

    return 0;
}

