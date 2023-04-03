#include "Helper.h"

uint64_t Eratosthenes_Basic(const uint64_t Nmax, bool vPrimes[], void*, void*)
{
    uint64_t i, j, step;
    uint64_t Nsqrt = (uint64_t)sqrt(Nmax);

    vPrimes[0] = vPrimes[1] = false;    // 0 and 1 are not primes
    for (i = 2; i <= Nmax; i++) vPrimes[i] = true;

    unsigned numPrimes = 1;             // accounting for 2
    AddPrime(2);

    for (i = 3; i <= Nsqrt; i += 2)
        if (vPrimes[i] == 1)
        {
            numPrimes++;                // count root primes
            AddPrime(i);

            for (j = i * i, step = 2 * i; j <= Nmax; j += step)
                vPrimes[j] = false; // flag non-primes
        }

    for ( /*continue counting primes*/; i <= Nmax; i += 2)
        if (vPrimes[i])
        {
            numPrimes++;
            AddPrime(i);
        }

    return numPrimes;
}

uint64_t Trial_Division(const uint64_t Nmax, uint64_t vPrimes[], void*, void*)
{
    uint64_t i, j;
    uint64_t Nsqr = (uint64_t)sqrt(Nmax);

    unsigned numPrimes = 1;             // accounting for 2
    AddPrime(2);

    for (i = 3; i <= Nmax; i += 2)
    {
        bool is_prime = true;
        uint64_t nsqrt = (uint64_t)ceil(sqrt(i));
        for (j = 1; (j < numPrimes) && (vPrimes[j] <= nsqrt); j++)
            if ((i % vPrimes[j]) == 0) { is_prime = false; break; }

        if (is_prime)
        {
            vPrimes[numPrimes++] = i;
            AddPrime(i);
        }
    }

    return numPrimes;
}

uint64_t Singleton_optimized311(uint64_t Nmax, uint64_t IP[], uint64_t IQ[], int32_t JQ[])
{
    uint64_t Nsqr = (uint64_t)ceil(sqrt(Nmax));

    IP[0] = 2; IP[1] = 3; IP[2] = 5; IP[3] = 7;
    uint64_t numPrimes = 4;             // accounting for 2, 3, 5, 7
    AddPrime(2); AddPrime(3); AddPrime(5); AddPrime(7);

    uint64_t iqi, n, step, nq, j, ij;
    int32_t jqi;
    iqi = 5 * 5; jqi = -2 * 5;          // initial values for prime number 5
    IQ[1] = 7 * 7; JQ[1] = 2 * 7;       // initial values for prime number 7
    nq = 1;                             // number of elements in the heap

    for (n = 11, step = 2; n <= Nmax; n += step, step = 6 - step)
    {
        if (n == iqi)
            do
            {
                //advance iqi
                iqi = (jqi > 0) ? (iqi + jqi + jqi) : (iqi - jqi);
                jqi = -jqi;
                if (iqi > IQ[1])
                {
                    //change iqi if necessary
                    std::swap(IQ[1], iqi); std::swap(JQ[1], jqi);
                    // push down the value if necessary
                    for (j = 2, ij = 1; (j <= nq); ij = j, j *= 2)
                    {
                        //j is the first child of ij; slide right if required
                        if ((j < nq) && (IQ[j] > IQ[j + 1])) j++;
                        if (IQ[j] < IQ[ij])
                        {
                            std::swap(IQ[j], IQ[ij]); std::swap(JQ[j], JQ[ij]);
                        }
                        else //we found the right slot
                            break;
                    }
                }
            } while (n == iqi);
        else
        {
            IP[numPrimes++] = n; AddPrime(n);

            if (n <= Nsqr) //if a root prime
            {
                nq++; //add the new values to the end of the heap 
                IQ[nq] = n * n;
                JQ[nq] = (uint32_t)(((n - (n / 3) * 3) == 1) ? 2 * n : -2 * n);
                // push up the value if necessary
                for (j = nq, ij = j / 2; (j > 1) && (IQ[j] < IQ[ij]); j = ij, ij /= 2)
                {
                    //ij is the father of j
                    std::swap(IQ[j], IQ[ij]); std::swap(JQ[j], JQ[ij]);
                }
                if (iqi > IQ[1])
                {
                    //change also iqi if necessary
                    std::swap(IQ[1], iqi); std::swap(JQ[1], jqi);
                }
            }
        }
    };

    return numPrimes;
}

uint64_t Mairson_Standard(uint64_t Nmax, uint64_t rlink[], uint64_t llink[], uint64_t DELETE[])
{
    auto crossoff = [&](uint64_t i)
    {
        rlink[llink[i]] = rlink[i];
        llink[rlink[i]] = llink[i];
    };

    uint64_t i;
    for (i = 1; i <= Nmax; i++)
    {
        rlink[i] = i + 1;
        llink[i] = i - 1;
        DELETE[i] = 0;
    }
    rlink[Nmax] = llink[1] = 0;

    uint64_t prime = 2;
    uint64_t nsqrt = (uint64_t)floor(sqrt(Nmax));
    while (prime <= nsqrt)
    {
        uint64_t factor = prime;
        uint64_t pointer = 0;
        uint64_t k = prime * factor;
        while (k <= Nmax)
        {
            pointer++;
            DELETE[pointer] = k;
            factor = rlink[factor];
            k = prime * factor;
        }
        for (i = 1; i <= pointer; i++)
            crossoff(DELETE[i]);
        prime = rlink[prime];
    }

    uint64_t numPrimes = 0;
    for (i = rlink[1]; i > 0; i = rlink[i])
    {
        AddPrime(i);
        numPrimes++;
    }
    return numPrimes;
}

uint64_t GalePratt(uint64_t Nmax, uint8_t sv[], void*, void*)
{
    const uint8_t BIT_MASK[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

    auto SetBit = [&](uint64_t bitidx)
    {
        bitidx /= 2;
        uint64_t idx = bitidx / 8;
        uint8_t bit = (uint8_t)(bitidx % 8);
        sv[idx] |= BIT_MASK[bit];
    };
    auto GetBit = [&](uint64_t bitidx)
    {
        bitidx /= 2;
        uint64_t idx = bitidx / 8;
        uint8_t bit = (uint8_t)(bitidx % 8);
        return (sv[idx] & BIT_MASK[bit]);
    };

    std::vector<uint64_t> s = { 1 };

    uint64_t i, k, m;
    for (i = 0; i <= Nmax; i++) sv[i] = 0;     // assume all primes

    for (i = 3; i <= Nmax / 2; i += 2)
        if (!GetBit(i))
        {
            std::vector<uint64_t> ts;
            for (auto j : s)
            {
                k = j;
                while (k < Nmax)
                {
                    if ((k != 1) and (k != i) and (k != j))
                        SetBit(k);
                    m = k * i;
                    if (m < Nmax)
                        ts.push_back(k);
                    k = m;
                }
            }
            s = ts;
        }

    uint64_t numPrimes = 1; AddPrime(2);
    for (i = 3; i <= Nmax; i += 2)
        if (!GetBit(i))
        {
            numPrimes++;
            AddPrime(i);
        }

    return numPrimes;
}

uint64_t GriesMisra_Standard(const uint64_t Nmax, uint64_t rlink[], uint64_t llink[], void*)
{
    auto crossoff = [&](uint64_t i)
    {
        rlink[llink[i]] = rlink[i];
        llink[rlink[i]] = llink[i];
    };

    uint64_t i;
    for (i = 1; i <= Nmax; i++)
    {
        rlink[i] = i + 1;
        llink[i] = i - 1;
    }
    rlink[Nmax] = llink[1] = 0;

    uint64_t nsqrt = (uint64_t)ceil(sqrt(Nmax));
    uint64_t p = 2;
    while (p <= nsqrt)
    {
        uint64_t q = p;
        uint64_t pq = p * q;
        while (pq <= Nmax)
        {
            uint64_t x = pq;
            while (x <= Nmax)
            {
                crossoff(x);
                x *= p;
            }
            q = rlink[q];
            pq = p * q;
        }
        p = rlink[p];
    }

    unsigned numPrimes = 0;
    for (i = rlink[1]; i > 0; i = rlink[i])
    {
        AddPrime(i);
        numPrimes++;
    }
    return numPrimes;
}

uint64_t SoP(const uint64_t N, uint64_t s[], void*, void*)

{
    for (uint64_t i = 0; i <= N; i++) s[i] = 0;

    uint64_t maxS = 1;
    uint64_t length = 2;
    uint64_t p = 3;

    auto next = [&](uint64_t w) { return s[w]; };
    auto prev = [&](uint64_t w) { return s[w - 1]; };

    auto Append = [&](uint64_t w)
    {
        s[maxS] = w;
        s[w - 1] = maxS;
        maxS = w;
    };

    auto Delete = [&](uint64_t pf)
    {
        uint64_t temp1 = s[pf - 1];
        uint64_t temp2 = s[pf];
        s[temp1] = temp2;
        s[temp2 - 1] = temp1;

    };

    auto ExtendTo = [&](uint64_t n)
    {
        uint64_t w = 1;
        uint64_t x = length + 1;

        while (x <= n)
        {
            Append(x);
            w = next(w);
            x = length + w;
        }
        length = n;
        if (length == N)
            Append(N + 2);
    };

    auto DeleteMultiples = [&](uint64_t p)
    {
        uint64_t f = p;

        while (p * f <= length)
            f = next(f);
        while (f > 1)
        {
            f = prev(f);
            Delete(p * f);
        }
    };

    while (p * p <= N)
    {
        if (length < N)
            ExtendTo(std::min(p * length, N));
        DeleteMultiples(p);
        p = next(1);
    }
    if (length < N)
        ExtendTo(N);

    //get the primes
    cChecker ckr;
    assert(ckr.check_next_prime(2));
    uint64_t  numPrimes = 1;    //account for 2
    AddPrime(2);
    p = 3;
    while (p <= N)
    {
        numPrimes++; AddPrime(p);
        p = next(p);
    }
    return numPrimes;
}

uint64_t Atkin(const uint64_t limit, bool sieve[], void*, void*)
{
    for (tpPrime i = 0; i < limit; i++) sieve[i] = false;

    // 2 and 3 are not generated by the alghorithm
    sieve[2] = true; sieve[3] = true;

    for (tpPrime x = 1; x * x < limit; x++)
    {
        for (tpPrime y = 1; y * y < limit; y++)
        {	// main loop - remove quadratic forms 
            tpPrime n = (4 * x * x) + (y * y);
            if (n <= limit && (n % 12 == 1 || n % 12 == 5))
                sieve[n] ^= true;

            n = (3 * x * x) + (y * y);
            if (n <= limit && n % 12 == 7)
                sieve[n] ^= true;

            n = (3 * x * x) - (y * y);
            if (x > y && n <= limit && n % 12 == 11)
                sieve[n] ^= true;
        }
    }
    // secondary loop - remove multiple of squares
    for (tpPrime i = 5; i * i < limit; i++) {
        if (sieve[i])
        {
            for (tpPrime j = i * i; j < limit; j += i * i)
                sieve[j] = false;
        }
    }
    tpPrime numprimes = 0;
    for (tpPrime i = 1; i < limit; i++)
        if (sieve[i])
        {
            numprimes++;
            AddPrime(i);
        }

    return numprimes;
}

constexpr uint64_t szM = 500'000;
bool vM[szM+1];

uint64_t Sieve357(const uint64_t Nmax, uint64_t IP[], uint64_t IQ[], int32_t JQ[], bool initialize = false)
{
    uint64_t numPrimes = 0;     //for this call
    for (int i = 0; i < szM; i++) vM[i] = true;

    static uint64_t nstart, nend;
    static uint32_t nq, stepq;
    static int32_t iq;      //idx of the last root prime generated so far
    if (initialize)
    {
        vM[1] = false;
        JQ[0] = 3;  IQ[0] = 9;
        iq = 0; stepq = 2; nq = 5;
        nstart = 0;
        nend = szM;
        numPrimes = 1; AddPrime(2); // account for 2
    }

    uint64_t nsqrt = (uint64_t)ceil(sqrt(nstart + szM));
    for (; nq <= nsqrt; nq += stepq, stepq = 6 - stepq)
    {   //generate rest of root primes for this iteration
        bool isprime = true;
        for (int i = 2; i <= iq; i++)
        {   //trial-division
            if ((nq / JQ[i]) * JQ[i] == nq)
            {
                isprime = false;
                break;
            }
        }
        if (isprime)
        {
            iq++;
            JQ[iq] = nq;
            IQ[iq] = ((uint64_t)nq) * nq;
        }
    }

    //strike out all composites
    for (int i = 0; i <= iq; i++)
    {
        uint64_t n = IQ[i];
        uint32_t stp = 2 * JQ[i];
        while (n < nend)
        {
            vM[n - nstart] = false;
            n += stp;
        }
        IQ[i] = n;
    }
    //move primes in IP
    for (int i = 1; i < szM; i += 2)
        if ((i + nstart) > Nmax)
            break;
        else
        if (vM[i])
            AddPrime(IP[numPrimes++] = nstart + i);

    nstart = nend;
    nend += szM;

    return numPrimes;
}

void Try357(const uint64_t LIMIT)
{
    nln();
    uint64_t sz = PI_Nmax(szM);
    uint64_t* vulPrimes = new uint64_t[sz + 1];
    //std::cout << sz << " : ";
    sz = PI_Nmax((uint64_t)ceil(sqrt(std::max(LIMIT, szM))));
    //std::cout << sz;
    uint64_t* IQ = new uint64_t[sz + 1];
    int32_t* JQ = new int32_t[sz + 1];

    nln(true);
    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();
    std::cout << " - Sieve 357 - ";
    for (auto i : Range<5>())
    {
        nln();

        uint64_t numPrimes = Sieve357(LIMIT, vulPrimes, IQ, JQ, true);
        for (uint64_t start = szM; start < LIMIT; start += szM)
            numPrimes += Sieve357(LIMIT, vulPrimes, IQ, JQ);

        std::cout << numPrimes << " primes up to " << LIMIT;
        times.push_back(tmr.LapTime(true));
    }
    nln(true);
    tmr.Stop();
    std::cout << "Average compute time: " << Average(times);

    nln(true);
    times.clear();

    delete[] vulPrimes;
    delete[] IQ;
    delete[] JQ;
}

uint64_t* p;
uint64_t* lpf;

tpPrime Ben(uint64_t Nmax)
{
    p[0] = lpf[0] = 0xFFFF; //idx 0 not used
    lpf[1] = lpf[2] = 0;
    uint64_t sz_lpf = 2;
    uint64_t sz_p = 0;

    for (uint64_t n = 2; n <= Nmax; n++)
    {
        if (lpf[n] == 0)
        {
            AddPrime(p[++sz_p] = n);
            lpf[n] = sz_p;
        }
        else
        {
            uint64_t q = n / p[lpf[n]];
            if (lpf[n] < lpf[q])
            {
                uint64_t r = lpf[n] + 1;
                lpf[q * p[r]] = r;
            }
        }
        lpf[++sz_lpf] = 0; lpf[++sz_lpf] = 1;
    }

    return sz_p;
}

void TryBen(const uint64_t LIMIT)
{
    uint64_t sz = PI_Nmax(LIMIT);
    p = new uint64_t[sz];
    lpf = new uint64_t[2 * LIMIT + 1];

    nln(true);
    cTimer tmr;
    std::vector<decltype(tmr.LapTime())> times;

    tmr.Start();
    std::cout << "Bengelloun sieve\n";
    for (auto i : Range<5>())
    {
        nln();

        auto numPrimes = Ben(LIMIT);

        std::cout << numPrimes << " primes up to " << LIMIT;
        times.push_back(tmr.LapTime(true));
    }
    nln(true);
    tmr.Stop();
    std::cout << "Average compute time: " << Average(times);
}
