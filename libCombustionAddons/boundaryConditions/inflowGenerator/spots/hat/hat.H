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

SourceFiles
    hatI.H
    hat.C
    hatIO.C

\*---------------------------------------------------------------------------*/

#ifndef hat_H
#define hat_H

#include "scalar.H"
#include "vector.H"
#include "Random.H"
#include "IOdictionary.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class hat;

Ostream& operator<<(Ostream&, const hat&);
Istream& operator>>(Istream&, hat&);

/*---------------------------------------------------------------------------*\
                           Class hat Declaration
\*---------------------------------------------------------------------------*/

class hat
{

public:

    class Parameters
    {
    public:
        scalar L_;    // integral length scale
        scalar Lspot_; // spot length

        Parameters();
        Parameters(const dictionary&);
        Parameters(scalar);

        void write(Ostream&) const;
    };

private:
    // Private data

    vector location_;
    vector epsilon_;

    // Private Member Functions

public:


    static scalar calcInfluenceLength(const Parameters& p);


    // Constructors
    hat();

    hat(Istream&);

    //- Construct from components
    hat(const vector&);
    
    //- Construct as copy
    hat(const hat&);
    
    // Selectors

    //- Select null constructed
    static autoPtr<hat> New(Istream&);


    // Destructor

    ~hat();
    

    // Member Functions

    // Access
    inline void setLocation(const vector& x) { location_=x; }
    inline const vector& location() const { return location_; }

    vector fluctuation(const Parameters&, const vector&) const;
    
    // Check

    // Edit
    void randomize(Random&);
    void moveForward(vector);
    
    // Write


    // Member Operators

    void operator=(const hat&);
    bool operator!=(const hat&) const;

    autoPtr<hat> clone() const;

    // Friend Functions

    // Friend Operators
    friend Ostream& operator<<(Ostream&, const hat&);
    friend Istream& operator>>(Istream&, hat&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
