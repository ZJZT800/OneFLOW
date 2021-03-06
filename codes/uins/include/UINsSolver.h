/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#pragma once
#include "HXDefine.h"
#include "INsSolver.h"
#include <map>
using namespace std;
BeginNameSpace( ONEFLOW )

class InvFlux;

class UINsSolver : public INsSolver
{
public:
    UINsSolver();
    ~UINsSolver();
public:
    IMPLEMENT_SOLVER_CLONE( UINsSolver )
public:
    void Run();
    void Init();
//public:
//    void ConvDiscret();
//    void DiffuDiscret();
public:
    void StaticInit();
//public:
    //void uINsSolver();
//public:
//    RealField* q;
//    MRField* rho;
//    MRField* mu;
};

//void INsPreflux();
//void INsCmpInv(); 
//void INsCmpVis();
//void INsCmpUnstead(); 
//void INsCmpSrc(); 
//void INsMomPre(); 
//void INsCmpFaceflux(); 
//void INsCorrectPresscoef(); 
//void INsCmpPressCorrectEquandUpdatePress();  
//void INsCmpSpeedCorrectandUpdateSpeed();  
//void INsUpdateFaceflux();   
//void INsUpdateRes();

EndNameSpace