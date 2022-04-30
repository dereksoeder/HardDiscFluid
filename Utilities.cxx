#include <tuple>
#include "Utilities.hxx"


std::mt19937_64 g_PRNG{std::random_device{}()};

std::mt19937_64 & GetPRNG()
{
    return g_PRNG;
}

double RandomReal()
{
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(g_PRNG);
}

double MaxwellBoltzmannSpeed(double M, double T)
{
#if 0
    // Maxwell-Boltzmann speed distribution is f(v) = (M/T) v exp(-(M*v*v/2) / T)
    // the peak occurs at v_peak = sqrt(T/M) and is f(sqrt(T/M)) = sqrt(M / T / e) ~ sqrt(M/T) * 0.60653066
    // unnormalized, f'(v) = v exp(-(M*v*v/2) / T) and its peak value is sqrt(T / M / e)

    double vofpeak = std::sqrt(T/M);
    double fatpeak = vofpeak * 0.60653066;  // note: this is for the unnormalized function f'(v)

    double Mover2T = (M/2.) / T;

    double s;
    do
    {
        s = RandomReal() * 4. * vofpeak;  // covers > 99.96% of distribution
    }
    while (RandomReal() * fatpeak > s * std::exp(-(s*s) * Mover2T));

    return s;
#else
    auto vxvy = MaxwellBoltzmannVelocity(M, T);
    return std::sqrt( std::get<0>(vxvy)*std::get<0>(vxvy) + std::get<1>(vxvy)*std::get<1>(vxvy) );
#endif
}

std::pair<double,double> MaxwellBoltzmannVelocity(double M, double T)
{
    std::normal_distribution<double> dist(0., std::sqrt(T/M));
    return std::make_pair(dist(GetPRNG()), dist(GetPRNG()));
}
