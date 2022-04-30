#pragma once

#include <memory>
#include <set>
#include <stdexcept>
#include <vector>
#include "Disc.hxx"


typedef enum
{
    LeftBoundary = 0,
    RightBoundary,
    TopBoundary,
    BottomBoundary,

    NumBoundaryPositions
} BoundaryPosition;

typedef enum
{
    OpenBoundary = 0,
    ReflectiveBoundary,
    PeriodicBoundary,
    ErraticBoundary,
    HotBoundary,

    NumBoundaryTypes
} BoundaryType;


#define BoundaryFlag_None 0U
#define LeftFlag   0x1U
#define RightFlag  0x2U
#define TopFlag    0x4U
#define BottomFlag 0x8U

class Arena
{
public:
    Arena(double width, double height,
          BoundaryType leftBC, BoundaryType rightBC, BoundaryType topBC, BoundaryType bottomBC)
      : m_Width(width),
        m_Height(height)
    {
        if (width <= 0.)
            throw std::invalid_argument{"width must be positive"};

        if (height <= 0.)
            throw std::invalid_argument{"height must be positive"};

        if ((leftBC < 0) || (leftBC >= NumBoundaryTypes))
            throw std::invalid_argument{"invalid left boundary type"};

        if ((rightBC < 0) || (rightBC >= NumBoundaryTypes))
            throw std::invalid_argument{"invalid right boundary type"};

        if ((topBC < 0) || (topBC >= NumBoundaryTypes))
            throw std::invalid_argument{"invalid top boundary type"};

        if ((bottomBC < 0) || (bottomBC >= NumBoundaryTypes))
            throw std::invalid_argument{"invalid bottom boundary type"};

        if (((leftBC == PeriodicBoundary) != (rightBC == PeriodicBoundary)) || ((topBC == PeriodicBoundary) != (bottomBC == PeriodicBoundary)))
            throw std::invalid_argument{"periodic boundaries must be paired"};

        m_BoundaryConditions[LeftBoundary]   = leftBC;
        m_BoundaryConditions[RightBoundary]  = rightBC;
        m_BoundaryConditions[TopBoundary]    = topBC;
        m_BoundaryConditions[BottomBoundary] = bottomBC;
    }

    void SetHotBoundaryTemperature(BoundaryPosition which, double temp)
    {
        if ((which < 0) || (which >= NumBoundaryPositions))
            throw std::invalid_argument{"invalid boundary position when setting hot boundary tempoerature"};

        if (m_BoundaryConditions[which] != HotBoundary)
            throw std::invalid_argument{"cannot set temperature for non-hot boundary"};

        m_HotBoundaryTemperatures[which] = temp;
    }

    double GetWidth() const
    { return m_Width; }

    double GetHeight() const
    { return m_Height; }

    double GetCurrentTime() const
    { return m_CurrentTime; }

    std::vector<Disc>::const_iterator GetDiscs() const
    { return m_Discs.cbegin(); }

    std::vector<Disc>::const_iterator GetDiscsEnd() const
    { return m_Discs.cend(); }

    std::vector<Disc>::size_type GetNumDiscs() const
    { return m_Discs.size(); }

    double GetInflowVx() const
    { return m_Inflow_Vx; }

    double GetInflowT() const
    { return m_Inflow_T; }

    void SetInflow(double Vx, double T)
    {
        m_Inflow_Vx = Vx;
        m_Inflow_T  = T;
    }

    inline double GetTotalWallMomentum() const
    { return m_CumulativeWallMomentum; }

    bool CanPlaceDisc(double x, double y, double r) const
    {
        if (r <= 0.)
            throw std::invalid_argument{"disc radius must be positive"};

        if (x <= r)
        {
            switch (m_BoundaryConditions[LeftBoundary])
            {
            case ReflectiveBoundary:
            case ErraticBoundary:
            case HotBoundary:
                return false;
            case PeriodicBoundary:
                if (!this->CanPlaceDiscUnbounded(x + m_Width, y, r))
                    return false;
                break;
            }
        }

        if (x >= (m_Width - r))
        {
            switch (m_BoundaryConditions[RightBoundary])
            {
            case ReflectiveBoundary:
            case ErraticBoundary:
            case HotBoundary:
                return false;
            case PeriodicBoundary:
                if (!this->CanPlaceDiscUnbounded(x - m_Width, y, r))
                    return false;
                break;
            }
        }

        if (y <= r)
        {
            switch (m_BoundaryConditions[TopBoundary])
            {
            case ReflectiveBoundary:
            case ErraticBoundary:
            case HotBoundary:
                return false;
            case PeriodicBoundary:
                if (!this->CanPlaceDiscUnbounded(x, y + m_Height, r))
                    return false;
                break;
            }
        }

        if (y >= (m_Height - r))
        {
            switch (m_BoundaryConditions[BottomBoundary])
            {
            case ReflectiveBoundary:
            case ErraticBoundary:
            case HotBoundary:
                return false;
            case PeriodicBoundary:
                if (!this->CanPlaceDiscUnbounded(x, y - m_Height, r))
                    return false;
                break;
            }
        }

        return this->CanPlaceDiscUnbounded(x, y, r);
    }

    bool Place(const Disc & disc)
    {
        if (!this->CanPlaceDisc(disc.x, disc.y, disc.R))
            return false;

        m_Discs.push_back(disc);
        m_DiscCollisions.emplace_back();
        m_DiscWrapFlags.push_back(this->GetWrapFlags(disc.x, disc.y, disc.R));

        this->UpdateCollisionsDisc(m_Discs.size() - 1);
        this->UpdateNextCollision();

        return true;
    }

    bool Step(double dtMax = 0.)
    {
        if (!m_NextCollision.IsValid())
        {
            if (dtMax > 0.)
            {
                this->Advance(dtMax);
                return true;
            }

            // nothing is going to happen next; caller would need to step by a definite time interval
            return false;
        }

        auto tnextcoll = (m_NextCollision.AbsoluteTime - m_CurrentTime);

        if ((dtMax < tnextcoll) && (dtMax > 0.))
        {
            // end of requested time interval < next collision

            this->Advance(dtMax);
            return true;
        }

        // collision is next

        this->Advance(tnextcoll);

        if (m_NextCollision.Type1 != Object_Disc)
            throw std::runtime_error{"first object in two-body collision record is not a disc!"};

        auto & disc = m_Discs[m_NextCollision.Index1];

        switch (m_NextCollision.Type2)
        {
        case Object_Disc:
            {
                auto & disc2 = m_Discs[m_NextCollision.Index2];

                if (m_NextCollision.Wrap2 == BoundaryFlag_None)
                {
                    disc.Collide(disc2);
                }
                else
                {
                    auto xold = disc2.x, yold = disc2.y;
                    disc2.x += ((m_NextCollision.Wrap2 & LeftFlag) ? m_Width  : (m_NextCollision.Wrap2 & RightFlag)  ? -m_Width  : 0.);
                    disc2.y += ((m_NextCollision.Wrap2 & TopFlag)  ? m_Height : (m_NextCollision.Wrap2 & BottomFlag) ? -m_Height : 0.);

                    disc.Collide(disc2);

                    disc2.x = xold;
                    disc2.y = yold;
                }
            }
            break;

        case Object_Boundary:
            if (m_NextCollision.Index2 >= NumBoundaryPositions)
                throw std::runtime_error{"disc is colliding with an invalid boundary!"};

            switch (m_BoundaryConditions[m_NextCollision.Index2])
            {
            case ReflectiveBoundary:
                {
                    if (!std::isfinite(disc.M))
                        throw std::runtime_error{"disc with infinite mass is colliding with wall!"};

                    auto horizline = (m_NextCollision.Index2 == TopBoundary) || (m_NextCollision.Index2 == BottomBoundary);
                    double vold = (horizline ? disc.vy : disc.vx);

                    disc.CollideBoundary(horizline);

                    m_CumulativeWallMomentum += std::abs((vold - (horizline ? disc.vy : disc.vx)) * disc.M);
                }
                break;

            case PeriodicBoundary:
                switch (m_NextCollision.Index2)
                {
                case LeftBoundary:
                    if (disc.x >= 0.)      m_DiscWrapFlags[m_NextCollision.Index1] |= LeftFlag;
                    break;
                case RightBoundary:
                    if (disc.x < m_Width)  m_DiscWrapFlags[m_NextCollision.Index1] |= RightFlag;
                    break;
                case TopBoundary:
                    if (disc.y >= 0.)      m_DiscWrapFlags[m_NextCollision.Index1] |= TopFlag;
                    break;
                case BottomBoundary:
                    if (disc.y < m_Height) m_DiscWrapFlags[m_NextCollision.Index1] |= BottomFlag;
                    break;
                    }
                break;

            case ErraticBoundary:
                {
                    if (!std::isfinite(disc.M))
                        throw std::runtime_error{"disc with infinite mass is colliding with wall!"};

                    auto vxold = disc.vx, vyold = disc.vy;
                    auto vold = std::sqrt(vxold*vxold + vyold*vyold);

                    auto theta = RandomReal() * 2.*PI;
                    disc.vx = vold * std::cos(theta);
                    disc.vy = vold * std::sin(theta);

                    switch (m_NextCollision.Index2)
                    {
                    case LeftBoundary:   disc.vx =  std::abs(disc.vx); break;
                    case RightBoundary:  disc.vx = -std::abs(disc.vx); break;
                    case TopBoundary:    disc.vy =  std::abs(disc.vy); break;
                    case BottomBoundary: disc.vy = -std::abs(disc.vy); break;
                    }

                    m_CumulativeWallMomentum += disc.M * std::sqrt( ((disc.vx - vxold) * (disc.vx - vxold)) + ((disc.vy - vyold) * (disc.vy - vyold)) );  // not well-founded, but there's no clear prescription here
                }
                break;

            case HotBoundary:
                {
                    if (!std::isfinite(disc.M))
                        throw std::runtime_error{"disc with infinite mass is colliding with wall!"};

                    // "half-Maxwellian" distribution (Ishiwata et al., Eq. 3) with explicit mass would read:
                    //   f(vx,vy) = 2/sqrt(pi * (T/(M/2))**3) * vy * exp(-(M/2)(vx*vx + vy*vy)/T)
                    // its peak is f(0, sqrt(T/M)) = M / T sqrt(2 pi e) ~ 0.241971 M / T

                    auto vxold = disc.vx, vyold = disc.vy;
                    double vpara, vperp;

                    auto hbtemp = m_HotBoundaryTemperatures[m_NextCollision.Index2];
                    auto vvar   = hbtemp / (disc.M / 2.);
                    auto vwidth = std::sqrt(vvar);

                    do
                    {
                        vpara = ((RandomReal() * 2.) - 1.) * 4. * vwidth;  // this should cover 99.96% of distribution
                        vperp =   RandomReal()             * 4. * vwidth;
                    }
                    while ( RandomReal() * 0.24197072451914337 >  // normalize to peak of "half-Maxwellian" distribution, Eq. 3 of  Ishiwata et al. (peak is f(0,sqrt(t)) = 1/(sqrt(2 pi e) t))
                            vperp * std::exp(-(vpara*vpara + vperp*vperp) / vvar) / std::sqrt(PI * vvar) );

                    switch (m_NextCollision.Index2)
                    {
                    case LeftBoundary:   disc.vx =  vperp;               disc.vy =  vpara; break;
                    case RightBoundary:  disc.vx = -vperp;               disc.vy =  vpara; break;
                    case TopBoundary:    disc.vx =  vpara + m_Inflow_Vx; disc.vy =  vperp; break;
                    case BottomBoundary: disc.vx =  vpara + m_Inflow_Vx; disc.vy = -vperp; break;
                    }

                    m_CumulativeWallMomentum += disc.M * std::sqrt( ((disc.vx - vxold) * (disc.vx - vxold)) + ((disc.vy - vyold) * (disc.vy - vyold)) );  // not well-founded, but there's no clear prescription here
                }
                break;
            } //switch(b.c.)
            break;

        case Object_BoundaryCrossing:
            if ((m_NextCollision.Index2 >= NumBoundaryPositions) || (m_BoundaryConditions[m_NextCollision.Index2] != PeriodicBoundary))
                throw std::runtime_error{"unexpected or invalid boundary-crossing event!"};

            switch (m_NextCollision.Index2)
            {
            case LeftBoundary:   disc.x += m_Width;  m_DiscWrapFlags[m_NextCollision.Index1] ^= (LeftFlag | RightFlag);  break;  // one flag had better be set; assuming so, this will switch it to the other side
            case TopBoundary:    disc.y += m_Height; m_DiscWrapFlags[m_NextCollision.Index1] ^= (TopFlag  | BottomFlag); break;
            case BottomBoundary: disc.y -= m_Height; m_DiscWrapFlags[m_NextCollision.Index1] ^= (TopFlag  | BottomFlag); break;

            case RightBoundary:
                {
                    disc.x -= m_Width;
                    m_DiscWrapFlags[m_NextCollision.Index1] ^= (LeftFlag | RightFlag);

                    if ((m_Inflow_Vx >= 1.E-6) && (m_Inflow_T > 0.))
                    {
                        std::pair<double,double> vxvy;
                        do
                        {
                            vxvy = MaxwellBoltzmannVelocity(disc.M, m_Inflow_T);
                        }
                        while (std::get<0>(vxvy) <= 0.);

                        disc.vx = std::get<0>(vxvy);
                        disc.vy = std::get<1>(vxvy);

                        for (int tries = 0; tries < 20; tries++)  // don't see this in Ishiwata et al., but if I don't do this then the initial rarefaction of the wake loops around
                        {
                            auto ytry = 1.001*disc.R + (RandomReal() * (m_Height - 2.002*disc.R));
                            if (CanPlaceDisc(disc.x, ytry, disc.R))
                            {
                                disc.y = ytry;
                                break;
                            }
                        }
                    } //if(inflow)
                }
                break;
            }
            break;

        case Object_BoundaryDetach:
            if ((m_NextCollision.Index2 >= NumBoundaryPositions) || (m_BoundaryConditions[m_NextCollision.Index2] != PeriodicBoundary))
                throw std::runtime_error{"unexpected or invalid boundary-detach event!"};

            m_DiscWrapFlags[m_NextCollision.Index1] &= ~m_FlagsForBoundary[m_NextCollision.Index2];
            break;

        default:
            throw std::runtime_error{"disc is colliding with an invalid object!"};
        } //switch

        this->UpdateCollisionsDisc(m_NextCollision.Index1);

        switch (m_NextCollision.Type2)
        {
        case Object_Disc:
            this->UpdateCollisionsDisc(m_NextCollision.Index2);
            break;
        }

        this->UpdateNextCollision();

        return true;
    }

    void RecomputeBoundaryCollisions()
    {
        for (size_t idx = 0; idx < m_Discs.size(); idx++)
        {
            auto & coll = m_DiscCollisions[idx];

            for (auto it = coll.begin(); it != coll.end(); )
            {
                if ((it->Type == Object_Boundary) || (it->Type == Object_BoundaryCrossing) || (it->Type == Object_BoundaryDetach))
                    it = coll.erase(it);
                else it++;
            }

            RecordBoundaryCollisionsDisc(idx);
        }
    }

private:
    enum ObjectType : unsigned int  // requires at least 3 bits
    {
        Object_None = 0,
        Object_Disc,
        Object_Boundary,
        Object_BoundaryCrossing,
        Object_BoundaryDetach
    } ;

    /* C++ doesn't seem to like this?
    enum BoundaryFlags : unsigned int
    {
        BoundaryFlag_None = 0U,
        LeftFlag   = (1U << LeftBoundary),
        RightFlag  = (1U << RightBoundary),
        TopFlag    = (1U << TopBoundary),
        BottomFlag = (1U << BottomBoundary)
    } ;*/
    typedef unsigned int BoundaryFlags;  // see also `#define`s at top

    static constexpr const BoundaryFlags m_FlagsForBoundary[NumBoundaryPositions] = { LeftFlag, RightFlag, TopFlag, BottomFlag };

    struct CollisionRecord
    {
        double AbsoluteTime;
        unsigned int Index : 24;
        ObjectType Type : 3;
        BoundaryFlags Wrap : 4;
        bool OurWrap : 1;

        CollisionRecord()
          : AbsoluteTime(-INF), Index(0), Type(Object_None), Wrap(BoundaryFlag_None), OurWrap(false)
        { }

        CollisionRecord(double absoluteTime, unsigned int index, ObjectType type, BoundaryFlags wrap = BoundaryFlag_None, bool ourwrap = false)
          : AbsoluteTime(absoluteTime), Index(index), Type(type), Wrap(wrap), OurWrap(ourwrap)
        { }

        inline bool IsValid() const
        { return (this->Type != Object_None); }

        inline void Invalidate()
        {
            this->Type = Object_None;
        }

        inline bool operator<(const CollisionRecord & other) const
        {
            // apparently `std::multiset` also uses `<` to determine how many elements are "equal" when removing, so we do have to consider all fields
            return (this->AbsoluteTime < other.AbsoluteTime) ||
                   ( (this->AbsoluteTime == other.AbsoluteTime) &&
                     ( (this->Type < other.Type) ||
                       ( (this->Type == other.Type) &&
                         ( (this->Index < other.Index) ||
                           ( (this->Index == other.Index) &&
                             ( (this->Wrap < other.Wrap) ||
                               ( (this->Wrap == other.Wrap) &&
                                 !(this->OurWrap) && other.OurWrap
                               )
                             )
                           )
                         )
                       )
                     )
                   );
        }
    } ;

    struct TwoBodyRecord
    {
        double AbsoluteTime;
        unsigned int Index1 : 24;
        ObjectType Type1 : 8;
        unsigned int Index2 : 24;
        ObjectType Type2 : 3;
        BoundaryFlags Wrap2 : 5;   // which doppelganger of a disc near one or more periodic boundaries is involved in the next collision

        TwoBodyRecord()
          : AbsoluteTime(-INF), Index1(0), Type1(Object_None), Index2(0), Type2(Object_None), Wrap2(BoundaryFlag_None)
        { }

        inline bool IsValid() const
        { return (this->Type1 != Object_None); }

        inline void Invalidate()
        {
            this->Type1 = Object_None;
        }
    } ;

    const double m_Width, m_Height;
    BoundaryType m_BoundaryConditions[NumBoundaryPositions];
    double m_HotBoundaryTemperatures[NumBoundaryPositions] = {};

    std::vector<Disc> m_Discs;
    std::vector<std::multiset<CollisionRecord>> m_DiscCollisions;  // important: `std::multiset` is ordered
    std::vector<BoundaryFlags> m_DiscWrapFlags;

    double m_CurrentTime = 0.;

    TwoBodyRecord m_NextCollision;

    double m_Inflow_Vx = 0.;
    double m_Inflow_T  = 0.;

    double m_CumulativeWallMomentum = 0.;

    inline BoundaryFlags GetWrapFlags(double x, double y, double r) const
    {
        return
            (((x < r)            && (m_BoundaryConditions[LeftBoundary]   == PeriodicBoundary)) ? LeftFlag   : BoundaryFlag_None) |
            (((x > m_Width  - r) && (m_BoundaryConditions[RightBoundary]  == PeriodicBoundary)) ? RightFlag  : BoundaryFlag_None) |
            (((y < r)            && (m_BoundaryConditions[TopBoundary]    == PeriodicBoundary)) ? TopFlag    : BoundaryFlag_None) |
            (((y > m_Height - r) && (m_BoundaryConditions[BottomBoundary] == PeriodicBoundary)) ? BottomFlag : BoundaryFlag_None);
    }

    bool CanPlaceDiscUnbounded(double x, double y, double r) const
    {
        for (const auto & disc : m_Discs)
        {
            auto dx = (disc.x - x), dy = (disc.y - y);
            auto sep = r + disc.R;
            if ((dx*dx + dy*dy) <= sep*sep)
                return false;
        }

        return true;
    }

    void Advance(double dt)
    {
        m_CurrentTime += dt;

        for (size_t idx = 0; idx < m_Discs.size(); idx++)
            m_Discs[idx].Advance(dt);
    }

    void UpdateCollisionsDisc(size_t index)
    {
        const auto & disc = m_Discs[index];
        auto & colls = m_DiscCollisions[index];

        for (const auto & coll: colls)
        {
            switch (coll.Type)
            {
            case Object_Disc:
                if (m_DiscCollisions[coll.Index].erase( CollisionRecord(coll.AbsoluteTime, index, Object_Disc, coll.Wrap, (coll.Wrap != BoundaryFlag_None && !coll.OurWrap)) ) != 1)
                    throw std::runtime_error{"failed to remove disc-disc collision record"};
                break;
            }
        }

        colls.clear();

        for (unsigned int iother = 0; iother < static_cast<unsigned int>(m_Discs.size()); iother++)  // downcast due to type of `Index` field (note: doesn't truncate to 24 bits though)
        {
            if (iother == index) continue;

            const auto & other = m_Discs[iother];
            auto flags = m_DiscWrapFlags[iother];

            if (flags == BoundaryFlag_None)
            {
                auto t = disc.TimeToCollision(other);
                if (t >= 0.)
                {
                    colls.emplace(m_CurrentTime + t, iother, Object_Disc);
                    m_DiscCollisions[iother].emplace(m_CurrentTime + t, index, Object_Disc);
                }
            }
            else
            {
                // check for collisions between this disc and other disc as well as its doppelgangers

                double tbest = -INF;
                BoundaryFlags wrapbest = BoundaryFlag_None;
                auto tmpdisc = other;

                for (int h = -1; h <= 1; h++)
                {
                    if ( (h != 0) && !(flags & ((h < 0) ? LeftFlag : RightFlag)) ) continue;

                    for (int v = -1; v <= 1; v++)
                    {
                        if ( (v != 0) && !(flags & ((v < 0) ? TopFlag : BottomFlag)) ) continue;

                        tmpdisc.x = other.x - (m_Width  * h);  // NOTE: minus because e.g. being near left boundary (`h = -1`) applies `LeftFlag`, but corresponding doppelganger is off the *right* side (`+ m_Width`, not `- m_Width`)
                        tmpdisc.y = other.y - (m_Height * v);

                        auto t = disc.TimeToCollision(tmpdisc);
                        if ((t >= 0.) && (t > tbest))
                        {
                            tbest = t;
                            wrapbest = ((h < 0) ? LeftFlag : (h == 0) ? BoundaryFlag_None : RightFlag) | ((v < 0) ? TopFlag : (v == 0) ? BoundaryFlag_None : BottomFlag);
                        }
                    }
                }

                if (tbest >= 0.)
                {
                    colls.emplace(m_CurrentTime + tbest, iother, Object_Disc, wrapbest, false);
                    m_DiscCollisions[iother].emplace(m_CurrentTime + tbest, index, Object_Disc, wrapbest, (wrapbest != BoundaryFlag_None));
                }
            } //if-else(flags)
        } //for(iother)

        if (m_DiscWrapFlags[index] != BoundaryFlag_None)
        {
            // now check for other discs' collisions with this disc's doppelgangers

            auto flags = m_DiscWrapFlags[index];
            auto tmpdisc = disc;

            for (unsigned int iother = 0; iother < static_cast<unsigned int>(m_Discs.size()); iother++)  // downcast due to type of `Index` field (note: doesn't truncate to 24 bits though)
            {
                if (iother == index) continue;

                const auto & other = m_Discs[iother];

                double tbest = -INF;
                BoundaryFlags wrapbest = BoundaryFlag_None;

                for (int h = -1; h <= 1; h++)
                {
                    if ( (h != 0) && !(flags & ((h < 0) ? LeftFlag : RightFlag)) ) continue;

                    for (int v = -1; v <= 1; v++)
                    {
                        if ( (v == 0) ? (h == 0) : !(flags & ((v < 0) ? TopFlag : BottomFlag)) )
                            continue;

                        tmpdisc.x = disc.x - (m_Width  * h);  // NOTE: minus because e.g. being near left boundary (`h = -1`) applies `LeftFlag`, but corresponding doppelganger is off the *right* side (`+ m_Width`, not `- m_Width`)
                        tmpdisc.y = disc.y - (m_Height * v);

                        auto t = other.TimeToCollision(tmpdisc);
                        if ((t >= 0.) && (t > tbest))
                        {
                            tbest = t;
                            wrapbest = ((h < 0) ? LeftFlag : (h == 0) ? BoundaryFlag_None : RightFlag) | ((v < 0) ? TopFlag : (v == 0) ? BoundaryFlag_None : BottomFlag);
                        }
                    }
                }

                if (tbest >= 0.)
                {
                    colls.emplace(m_CurrentTime + tbest, iother, Object_Disc, wrapbest, false);
                    m_DiscCollisions[iother].emplace(m_CurrentTime + tbest, index, Object_Disc, wrapbest, true);
                }
            } //for(iother)
        } //if(flags)

        RecordBoundaryCollisionsDisc(index);
    }

    void RecordBoundaryCollisionsDisc(size_t index)
    {
        const auto & disc = m_Discs[index];
        auto & colls = m_DiscCollisions[index];

#define CHECKBOUNDARYCOLLISION(__which, __coord, __cmp, __pos, __sign, __horizontal, __flag) \
if (disc.__coord  __cmp  __pos) \
{ \
    if ((m_BoundaryConditions[__which] == ReflectiveBoundary) || (m_BoundaryConditions[__which] == ErraticBoundary) || (m_BoundaryConditions[__which] == HotBoundary)) \
    { \
        auto t = disc.TimeToBoundaryCollision(__pos, __horizontal); \
        if (t >= 0.) \
        { \
            colls.emplace(m_CurrentTime + t, __which, Object_Boundary); \
        } \
    } \
    else if (m_BoundaryConditions[__which] == PeriodicBoundary) \
    { \
        if (m_DiscWrapFlags[index] & __flag)  /* record boundary crossing if disc is already "wrapping" (overlapping periodic boundary), otherwise record "collision" with boundary */ \
        { \
            auto t = disc.TimeToBoundaryCrossing(__pos, __horizontal); \
            if (t >= 0.) \
            { \
                colls.emplace(m_CurrentTime + t, __which, Object_BoundaryCrossing); \
            } \
            else \
            { \
                t = disc.TimeToBoundaryCrossing(__pos  __sign  disc.R, __horizontal); \
                if (t >= 0.) \
                { \
                    colls.emplace(m_CurrentTime + t, __which, Object_BoundaryDetach); \
                } \
            } \
        } \
        else \
        { \
            auto t = disc.TimeToBoundaryCollision(__pos, __horizontal); \
            if (t >= 0.) \
            { \
                colls.emplace(m_CurrentTime + t, __which, Object_Boundary); \
            } \
        } \
    } \
}

        CHECKBOUNDARYCOLLISION(LeftBoundary,   x, >=, 0.,       +, false, LeftFlag);
        CHECKBOUNDARYCOLLISION(RightBoundary,  x, <=, m_Width,  -, false, RightFlag);
        CHECKBOUNDARYCOLLISION(TopBoundary,    y, >=, 0.,       +, true,  TopFlag);
        CHECKBOUNDARYCOLLISION(BottomBoundary, y, <=, m_Height, -, true,  BottomFlag);

#undef CHECKBOUNDARYCOLLISION
    }

    void UpdateNextCollision()
    {
        m_NextCollision.Invalidate();

        // only check disc collision lists; infinite line and line segment lists are redundant

        for (auto idx = 0; idx < m_DiscCollisions.size(); idx++)
        {
            const auto & colllist = m_DiscCollisions[idx];

            if (colllist.size() > 0)
            {
                const auto & coll = *colllist.cbegin();

                if (!m_NextCollision.IsValid() || (coll.AbsoluteTime < m_NextCollision.AbsoluteTime))
                {
                    m_NextCollision.AbsoluteTime = coll.AbsoluteTime;
                    m_NextCollision.Index1 = (coll.OurWrap ? coll.Index  : idx);  // if `idx`th (not `coll.Index`th) disc's doppelganger is colliding, put it in the `Xxx2` fields since they include a wrap field
                    m_NextCollision.Type1  = (coll.OurWrap ? coll.Type   : Object_Disc);
                    m_NextCollision.Index2 = (coll.OurWrap ? idx         : coll.Index);
                    m_NextCollision.Type2  = (coll.OurWrap ? Object_Disc : coll.Type);
                    m_NextCollision.Wrap2  = coll.Wrap;
                }
            }
        }
    }
} ;
