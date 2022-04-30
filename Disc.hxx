#pragma once

#include <cmath>
#include <stdexcept>
#include "Utilities.hxx"


extern double g_CylinderMomentum;  // instantiated in main.cxx

class Disc
{
public:
    double x, y;
    double vx, vy;
    double M, R;

    Disc(double x_, double y_, double R_, double vx_ = 0., double vy_ = 0., double M_ = 0.)
      : x(x_), y(y_), R(R_), vx(vx_), vy(vy_), M((M_ > 0.) ? M_ : PI*R_*R_)
    {
        if (R_ <= 0.)
            throw std::invalid_argument{"disc radius must be positive"};
    }

    void Advance(double dt)
    {
        this->x += (this->vx * dt);
        this->y += (this->vy * dt);
    }

    double TimeToBoundaryCollision(double pos, bool horizontalLine) const
    {
        auto v = (horizontalLine ? this->vy : this->vx);
        if (v == 0.) return -INF;  // no collision if not moving, even if already touching

        auto r  = (horizontalLine ? this->y : this->x);
        auto dr = (pos - r);

        if (std::abs(dr) <= this->R)
            return ((dr * v) > 0.) ? 0. : -INF;  // if already touching, immediate collision only if disc is headed toward boundary (`dr * v` like dot-product, just to check sign; not a meaningful quantity)

        if (dr > 0.) dr -= this->R;
        else         dr += this->R;

        return (v == 0.) ? -INF : (dr / v);
    }

    double TimeToBoundaryCrossing(double pos, bool horizontalLine) const
    {
        auto r  = (horizontalLine ?  this->y  : this->x);
        auto v  = (horizontalLine ? this->vy : this->vx);
        auto dr = (pos - r);

        return ((dr * v) == 0.) ? -INF : (dr / v);  // future (past) crossing if headed toward (away from) boundary but not currently there (`dr * v` like dot-product, just to check sign; not a meaningful quantity)
    }

    double TimeToCollision(const Disc & disc) const
    {
        auto sep   = (disc.R + this->R);
        auto sepsq = sep*sep;

        auto dx   = (disc.x - this->x), dy = (disc.y - this->y);
        auto drsq = (dx*dx + dy*dy);

        auto dvx    = (disc.vx - this->vx), dvy = (disc.vy - this->vy);
        auto bover2 = (dx*dvx + dy*dvy);

        auto c = (drsq - sepsq);

        if (c <= 1.E-9)  // if already touching, immediate collision if also headed toward each other, otherwise no collision so that they can separate (need small tolerance because solutions computed below can be imprecise)
        {
            return (bover2 < 0.) ? 0. : -INF;  // `bover2` = dr.dv; want it to be *negative* for collision, so that separation is antiparallel to net velocity (i.e., they're approaching)
        }

        auto dvsq = (dvx*dvx + dvy*dvy);

        // { -(dx*dvx + dy*dvy) +/- sqrt[ (dx*dvx + dy*dvy)**2 - (dvx**2 + dvy**2)*(dx**2 + dy**2 - sep**2) ] } / (dvx**2 + dvy**2)

        auto bsqover4 = bover2*bover2;
        auto ac       = dvsq * c;

        if (bsqover4 < ac)
            return -INF;  // no (real) solution: will never collide

        auto sqrtterm = std::sqrt(bsqover4 - ac);
        return ((-bover2 >= sqrtterm) ? (-bover2 - sqrtterm) : (-bover2 + sqrtterm)) / dvsq;  // return smaller positive solution (or no particular negative solution if both solutions are negative)
    }

    void Collide(Disc & disc)
    {
        auto thisfinite = std::isfinite(this->M), discfinite = std::isfinite(disc.M);

        if (!thisfinite && discfinite)
        {
            disc.Collide(*this);
            return;
        }

        auto dx   = (disc.x - this->x), dy = (disc.y - this->y);
        auto drsq = (dx*dx + dy*dy);

        auto sep   = (disc.R + this->R);
        auto sepsq = sep*sep;

        //if (drsq > sepsq)
        //    throw std::runtime_error{"cannot collide non-contacting discs"};

        if (thisfinite == discfinite)
        {
            double thisM = (thisfinite ? this->M : 1.), discM = (discfinite ? disc.M : 1.);  // if both have infinite mass, treat masses as equal

            auto pi1x = thisM * this->vx, pi1y = thisM * this->vy;
            auto pi2x = discM *  disc.vx, pi2y = discM *  disc.vy;

            auto j  = 2. * (discM * (pi1x * dx + pi1y * dy) - thisM * (pi2x * dx + pi2y * dy)) / (thisM + discM) / drsq;
            auto jx = j * dx, jy = j * dy;

            this->vx -= (jx / thisM);
            this->vy -= (jy / thisM);

            disc.vx  += (jx / discM);
            disc.vy  += (jy / discM);
        }
        else  // `(thisfinite && !discfinite)` ensured by earlier check and return; must temporarily boost to disc's rest frame to suppress its infinite momentum
        {
#if 0
            auto pi1x = this->M * (this->vx - disc.vx), pi1y = this->M * (this->vy - disc.vy);
            auto j = 2. * (pi1x * dx + pi1y * dy) / this->M / drsq;

            this->vx += disc.vx - (j * dx);
            this->vy += disc.vy - (j * dy);
#else
            // infinite-mass disc is the cylinder, has erratic boundary conditions

            auto r = std::sqrt(drsq);
            auto v = std::sqrt((this->vx * this->vx) + (this->vy * this->vy));

            auto nx = -dx / r, ny = -dy / r;  // normal vector from cylinder (`disc`) to this disc
            auto tx =  ny,     ty = -nx;      // tangent vector (sign won't matter)

            auto theta = RandomReal() * PI;
            auto costheta = std::cos(theta), sintheta = std::sin(theta);
            auto newvx = v * ((costheta * tx) + (sintheta * nx));
            auto newvy = v * ((costheta * ty) + (sintheta * ny));

            auto dvx = (newvx - this->vx);
            auto dvy = (newvy - this->vy);
            g_CylinderMomentum -= dvx * this->M;

            this->vx = newvx;
            this->vy = newvy;
#endif
        }
    }

    void CollideBoundary(bool horizontalLine)
    {
        if (horizontalLine)
             this->vy = -this->vy;
        else this->vx = -this->vx;
    }
} ;
