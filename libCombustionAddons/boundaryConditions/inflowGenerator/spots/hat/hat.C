/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    hat

Description

Author

\*----------------------------------------------------------------------------*/

#include "hat.H"

namespace Foam
{

hat::Parameters::Parameters
(
)
    :
    L_(0.0),    // integral length scale
    Lspot_(calcInfluenceLength(*this))
{
}

hat::Parameters::Parameters
(
    const dictionary& dict
)
    :
    L_(readScalar(dict.lookup("L"))),    // integral length scale
    Lspot_(calcInfluenceLength(*this))
{
}

hat::Parameters::Parameters
(
    scalar L    // integral length scale

):
    L_(L),    // integral length scale
    Lspot_(calcInfluenceLength(*this))
{
}

void hat::Parameters::write
(
    Ostream& os
) const
{
    os.writeKeyword("L")
        << L_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar hat::calcInfluenceLength(const Parameters& p)
{
#warning Please check the factor!
    return (4./3.)*p.L_;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hat::hat()
: 
    location_(pTraits<vector>::zero),
    epsilon_(pTraits<vector>::zero)
{
}

hat::hat
(
    Istream& s
)
: 
    location_(s),
    epsilon_(s)
{
    Info<<location_<<endl;
}

hat::hat(const vector& loc)
:
    location_(loc),
    epsilon_(pTraits<vector>::zero)
{}


hat::hat(const hat& o)
:
    location_(o.location_),
    epsilon_(o.epsilon_)
{}

vector hat::fluctuation(const Parameters& p, const vector& x) const
{
    vector delta_x = x - location_;

    if 
        (
            (delta_x.x()  < p.Lspot_ / 2.0) &&
            (delta_x.y()  < p.Lspot_ / 2.0) &&
            (delta_x.z()  < p.Lspot_ / 2.0)
        )
    {
      vector f=
           (1.0 - 2.0*delta_x.x()  / p.Lspot_)
          *(1.0 - 2.0*delta_x.y()  / p.Lspot_)
          *(1.0 - 2.0*delta_x.z()  / p.Lspot_)
          * pTraits<vector>::one;

      return cmptMultiply(epsilon_, f);
    }
  else
    return pTraits<vector>::zero;
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<hat> hat::New(Istream& s)
{
    Info<<"reading"<<endl;
    return autoPtr<hat>(new hat(s));
}


void hat::randomize(Random& rand)
{
    rand.randomise(epsilon_);
    epsilon_ -= pTraits<vector>::one*0.5;
    epsilon_ *= 2.0;
}

void hat::moveForward(vector delta)
{
    location_+=delta;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hat::~hat()
{}


autoPtr<hat> hat::clone() const
{
    return autoPtr<hat>
        (
            new hat(*this)
        );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void hat::operator=(const hat& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("hat::operator=(const hat&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    location_=rhs.location_;
    epsilon_=rhs.epsilon_;
}

bool hat::operator!=(const hat& o) const
{
    return 
        (location_!=o.location_)
        ||
        (epsilon_!=o.epsilon_);
}

Ostream& operator<<(Ostream& s, const hat& ht)
{
    s<<ht.location_<<endl;
    s<<ht.epsilon_<<endl;
    return s;
}

Istream& operator>>(Istream& s, hat& ht)
{
    vector loc(s);
    vector eps(s);
    ht.location_=loc;
    ht.epsilon_=eps;
    return s;
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //

}
