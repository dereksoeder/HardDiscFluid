#include <stdio.h>
#include <string.h>
#include <cmath>
#include "Arena.hxx"
#include "Utilities.hxx"



////////////////////////////////////////////////////////////////

// these tests are generally based on Ishiwata et al., Int. J. Mod. Phys. C 15 (2004) 1413-1424, https://doi.org/10.1142/S0129183104006820

double g_CylinderMomentum;  // extern in Disc.hxx


////////////////////////////////
// TestPlanePoiseuille
////////////////////////////////
int TestPlanePoiseuille(int argc, const char * argv[])
{
    double density = -1.;
    if (argc != 1 || sscanf(argv[0], "%lg", &density) != 1 || !std::isfinite(density) || density < 0.2 || density >= 0.500001)
    {
        fprintf(stderr, "TestPlanePoiseuille: density must be between 0.2 and 0.5\n");
        return 1;
    }

    constexpr double Lx = 300., Ly = 100.;
    constexpr double R = 1.5;
    constexpr double M = 9.0;
    constexpr double T = 2.0;
    constexpr double g = 0.001;

    const int ndiscs = static_cast<int>(std::floor(density * Lx * Ly / M));

    Arena arena(Lx, Ly, PeriodicBoundary, PeriodicBoundary, HotBoundary, HotBoundary);
    arena.SetHotBoundaryTemperature(TopBoundary,    T);
    arena.SetHotBoundaryTemperature(BottomBoundary, T);

    for (int n = 0; n < ndiscs; n++)
    {
        double x, y;
        do
        {
            x = 1.001*R + (RandomReal() * (arena.GetWidth()  - 2.002*R));
            y = 1.001*R + (RandomReal() * (arena.GetHeight() - 2.002*R));
        }
        while (!arena.CanPlaceDisc(x, y, R));

        auto vxvy = MaxwellBoltzmannVelocity(M, T);

        if (!arena.Place(Disc(x, y, R, std::get<0>(vxvy), std::get<1>(vxvy), M)))
            fprintf(stderr, "WARNING: placing disc #%d at (%lg, %lg) failed\n", n, x, y);
    }

    constexpr int nlanes = 30;

    double sums[nlanes]  = { };
    int nsamples[nlanes] = { };

    double nextsampletime = 2000.;
    double lastacceltime  = 0.;

    while (arena.GetCurrentTime() < 15000.)
    {
        if (arena.GetCurrentTime() >= nextsampletime)
        {
            nextsampletime += 1.;

            for (auto it = arena.GetDiscs(); it != arena.GetDiscsEnd(); it++)
            {
                auto iy = static_cast<int>(std::floor( nlanes * (it->y - R) / (arena.GetHeight() - 2.*R) ));
                if ((iy >= 0) && (iy < nlanes))
                {
                    sums[iy] += it->vx;
                    nsamples[iy]++;
                }
            }
        }

        arena.Step();

        double newtime = arena.GetCurrentTime();
        if ((newtime - lastacceltime) > 1.E-6)
        {
            for (auto it = arena.GetDiscs(); it != arena.GetDiscsEnd(); it++)
            {
                auto pdisc = const_cast<Disc*>(&*it);
                pdisc->vx += g * (newtime - lastacceltime);
            }

            lastacceltime = newtime;

            arena.RecomputeBoundaryCollisions();  // uniformly accelerating all discs doesn't change inter-disc next collision times, but does change next boudnary collision times
        }
    }

    for (int n = 0; n < nlanes; n++)
        printf("%.10lg %.10lg\n", R + (n * (arena.GetHeight() - 2.*R) / nlanes), sums[n] / nsamples[n]);

    return 0;
} //TestPlanePoiseuille


////////////////////////////////
// TestCylinder
////////////////////////////////
int TestCylinder(int argc, const char * argv[])
{
    int re = -1;
    if ((argc != 1) || (sscanf(argv[0], "%d", &re) != 1) || (re <= 0))
    {
        fprintf(stderr, "TestCylinder: Reynolds number must be positive\n");
        return 1;
    }

    int ndiscs;
    double Vx, d, Lx, Ly;
    double xcyl;

    switch (re)  // Ishiwata et al., Table 2  (I had to guess xcyl)
    {
  //case   1:  Vx = 0.012;  d = 150.;  Lx = 1000.;  Ly =  750.;  ndiscs = 292931;  xcyl =  400.;  break;
    case   6:  Vx = 0.11;   d = 100.;  Lx =  800.;  Ly =  800.;  ndiscs = 252858;  xcyl =  300.;  break;  // Fig. 5
    case  11:  Vx = 0.2;    d = 100.;  Lx =  800.;  Ly =  800.;  ndiscs = 252858;  xcyl =  300.;  break;
    case  33:  Vx = 1.;     d =  60.;  Lx =  800.;  Ly =  800.;  ndiscs = 254869;  xcyl =  200.;  break;  // Fig. 6
    case  66:  Vx = 1.;     d = 120.;  Lx = 1000.;  Ly =  800.;  ndiscs = 315476;  xcyl =  400.;  break;
    case 106:  Vx = 1.;     d = 200.;  Lx = 1600.;  Ly = 1200.;  ndiscs = 566575;  xcyl =  460.;  break;  // Fig. 7
    default:
        throw std::invalid_argument{"unsupported Reynolds number"};
    }

    constexpr double R = 1.5;
    constexpr double M = 9.0;
    ndiscs /= 9;

    constexpr double T = 2.;

    Arena arena(Lx, Ly, PeriodicBoundary, PeriodicBoundary, HotBoundary, HotBoundary);
    arena.SetHotBoundaryTemperature(TopBoundary,    T);
    arena.SetHotBoundaryTemperature(BottomBoundary, T);
    arena.SetInflow(Vx, T);

    arena.Place(Disc(xcyl, arena.GetHeight()/2., d/2., 0., 0., INF));  // infinite-mass disc is cylinder, has special no-slip condition in Disc.hxx

    for (int n = 0; n < ndiscs; n++)
    {
        double x, y;
        do
        {
            x = 1.001*R + (RandomReal() * (arena.GetWidth()  - 2.002*R));
            y = 1.001*R + (RandomReal() * (arena.GetHeight() - 2.002*R));
        }
        while (!arena.CanPlaceDisc(x, y, R));

        auto vxvy = MaxwellBoltzmannVelocity(M, T);

        if (!arena.Place(Disc(x, y, R, std::get<0>(vxvy), std::get<1>(vxvy), M)))
            fprintf(stderr, "WARNING: placing disc #%d at (%lg, %lg) failed\n", n, x, y);
    }

    int nsamples[120][160]   = { };
    double flow[120][160][2] = { };

    const int ncellsx = static_cast<int>(std::floor(Lx / 20.));
    const int ncellsy = static_cast<int>(std::floor(Ly / 20.));

    double startedtime = -1.;
    double nextsampletime = 1000.;

    while (arena.GetCurrentTime() < 5000.)
    {
        if (arena.GetCurrentTime() >= nextsampletime)
        {
            if (startedtime <= 0.)
            {
                startedtime = arena.GetCurrentTime();
                g_CylinderMomentum = 0.;
            }

            nextsampletime += 1.;

            for (auto it = arena.GetDiscs(); it != arena.GetDiscsEnd(); it++)
            {
                auto ix = static_cast<int>(std::floor( ncellsx *  it->x      /  arena.GetWidth() ));
                auto iy = static_cast<int>(std::floor( ncellsy * (it->y - R) / (arena.GetHeight() - 2.*R) ));

                if ((ix >= 0) && (ix < ncellsx) && (iy >= 0) && (iy < ncellsy))
                {
                    flow[iy][ix][0] += it->vx;
                    flow[iy][ix][1] += it->vy;
                    nsamples[iy][ix]++;
                }
            }
        }

        arena.Step();
    }

    printf("# %.10lg\n", g_CylinderMomentum / (arena.GetCurrentTime() - startedtime));

    for (int iy = 0; iy < ncellsy; iy++)
    {
        for (int ix = 0; ix < ncellsx; ix++)
            printf("%.10lg %.10lg ", flow[iy][ix][0] / nsamples[iy][ix], flow[iy][ix][1] / nsamples[iy][ix]);
        printf("\n");
    }

    return 0;
} //TestCylinderAverageFlow



////////////////////////////////////////////////////////////////

typedef int (* TestFunction)(int, const char * *);

static const struct
{
    TestFunction Func;
    const char * Label;
} g_Tests[] =
{
    { TestPlanePoiseuille, "Plane Poiseuille flow: velocity profile  (arg: density)" },
    { TestCylinder, "Flow around a cylinder: time-averaged flow and drag force  (arg: Re)" }
};

constexpr size_t NumTests = (sizeof(g_Tests) / sizeof(g_Tests[0]));


////////////////////////////////
// main
////////////////////////////////
int main(int argc, const char * argv[])
{
    int ntest = -1;
    if ((argc < 2) || (sscanf(argv[1], "%d", &ntest) != 1) || (ntest < 0) || (ntest >= NumTests))
    {
        fprintf(stderr, "Usage: %s testnumber [arg]\n\n"
                        "Where:  testnumber  is one of the following:\n\n", argv[0]);

        for (int i = 0; i < NumTests; i++)
            fprintf(stderr, "%7d : %s\n", i, g_Tests[i].Label);

        fprintf(stderr, "\n");
        return 1;
    }

    return g_Tests[ntest].Func(argc - 2, argv + 2);
} //main
