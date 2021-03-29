/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 31 "/home/sv015/OF/TVG-X-4-S78/0/U.boundaryField.inlet.#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_b54757a7e66dc9a109601977ce18d8c5012ac1f4
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 49 "/home/sv015/OF/TVG-X-4-S78/0/U.boundaryField.inlet.#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>
		(
                    dict.parent().parent()
                );

                const fvMesh& mesh = refCast<const fvMesh>(d.db());

                const label id = mesh.boundary().findPatchID("inlet");
                const fvPatch& patch = mesh.boundary()[id];

                //vectorField U(mesh.boundary()[id].size(), vector(0, 0, 0));
                vectorField U(patch.size(), vector(0, 0, 0));


                const scalar pi = constant::mathematical::pi;

                const scalar S = 0.78; 
                const scalar a = 1.0;

                const scalar z_0 =  0.00004;
                const scalar C_1 =  0.09875;

                const scalar Vin = 1;
		//const scalar Vin = 0.2;
                for (int i=0; i<patch.size(); i++)
                {
                	const scalar x = patch.Cf()[i][0];
			const scalar y = patch.Cf()[i][1];
			const scalar z = patch.Cf()[i][2];

                        const scalar U_r = Vin*C_1*(log(1+(z/z_0)));
                        const scalar U_t = 2*S*a*U_r;
                        
                        if (x>=0 && y>=0)	// 1st Quadrant
                        {
                        	const scalar alpha =  pi + atan(y/x);
                        	//const scalar theta =  atan(U_t/U_r);
                        	U[i] = vector (((U_r)*cos(alpha)+(U_t)*cos(alpha-pi/2)),((U_r)*sin(alpha)+(U_t)*sin(alpha-pi/2)),0.0);
			}
			else if(x<=0 && y>=0)	// 2nd Quadrant
			{	
				const scalar alpha =  atan(y/x);
                        	//const scalar theta =  atan(U_t/U_r);
                        	U[i] = vector (((U_r)*cos(alpha)+(U_t)*cos(alpha-pi/2)),((U_r)*sin(alpha)+(U_t)*sin(alpha-pi/2)),0.0);
			}
			else if(x<=0 && y<=0)	// 3rd Quadrant
			{	
				const scalar alpha =  2*pi + atan(y/x);
                        	//const scalar theta =  atan(U_t/U_r);
                        	U[i] = vector (((U_r)*cos(alpha)+(U_t)*cos(3*pi/2+alpha)),((U_r)*sin(alpha)+(U_t)*sin(3*pi/2+alpha)),0.0);
			}
			else	// 4th Qaudrant
			{	
				const scalar alpha =  pi + atan(y/x);
                        	//const scalar theta =  atan(U_t/U_r);
                        	U[i] = vector (((U_r)*cos(alpha)+(U_t)*cos(alpha-pi/2)),((U_r)*sin(alpha)+(U_t)*sin(alpha-pi/2)),0.0);
			}


		}

                writeEntry(os, "", U);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

