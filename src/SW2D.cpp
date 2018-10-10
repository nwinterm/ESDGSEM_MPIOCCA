#include "SW2D.h"
//#include "Constants.hpp"





SW2D::SW2D(const int tcase, const dfloat postol)
{
    Testcase = tcase;
    PosPresTOL = postol;

}
//  void InitQ(const int Nelem,const int ngl,const int ngl2,const dfloat x[NoSpaceDofs],const dfloat y[NoSpaceDofs], dfloat q[NoDofs],const dfloat t,const dfloat b[NoSpaceDofs]){

void SW2D :: InitQ(const int IsMeshSplit,const MeshPartitioning MeshSplit,const int Nelem,const int ngl,const int ngl2,
                   const dfloat x[],const dfloat y[], dfloat q[],const dfloat t,
                   const dfloat * b, const dfloat g_const, dfloat h_0)
{

    int globalEleID;

    for(int ie=0; ie<Nelem; ++ie)
    {
        if (IsMeshSplit==1)
        {
            globalEleID  = MeshSplit.MyElementLocalToGlobal[ie];
        }
        else
        {
            globalEleID  = ie+1;
        }
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                //erst alle Punkte die zu j=0 gehoeren, dann j=1, usw.
                int qid = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2  + j*ngl+i;


                dfloat qNodal[Neq];

                InitQNodal(x[xid],y[xid],qNodal,t,b[xid],g_const,h_0);




                if (Testcase==3)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 1.0;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
                        qNodal[0] = 0.1;
                    }
                    if( globalEleID % NelemX == 0)
                    {
                        qNodal[0] = 0.1;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 1.0;
                    }


                }
                if (Testcase==9)
                {
                    int NelemY = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID < (NelemY/2)*NelemY+1)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID > (NelemY/2)*NelemY+1)
                    {
                        qNodal[0] = 5.0;
                    }


                }

                if (Testcase==20)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
                        qNodal[0] = 5.0;
                    }
                    if( globalEleID % NelemX == 0)
                    {
                        qNodal[0] = 5.0;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 10.0;
                    }

                }
                if (Testcase==21)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 2.5;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
                        qNodal[0] = 1.5;
                    }
                    if( globalEleID % NelemX == 0)
                    {
                        qNodal[0] = 1.5;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 2.5;
                    }

                }

                if (Testcase==29)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
//                    qNodal[0] = max(pow(10.0,-6),0.1-b[xid]);
                        qNodal[0] = 5.0;
                    }
                    if( globalEleID % NelemX == 0)
                    {
//                    qNodal[0] = max(pow(10.0,-6),0.1-b[xid]);
                        qNodal[0] = 5.0;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 10.0;
                    }


                }

                if (Testcase==30)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
                        qNodal[0] = 0.0;
                    }
                    if( globalEleID % NelemX == 0)
                    {
                        qNodal[0] = 0.0;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 10.0;
                    }


                }
                if (Testcase==28)
                {
                    int NelemX = sqrt(MeshSplit.global_NumElements);
                    if( globalEleID % (NelemX/2) == 0)
                    {
                        qNodal[0] = 1.0;
                    }
                    if( globalEleID % (NelemX/2) == 1)
                    {
                        qNodal[0] = 0.0;
                    }
                    if( globalEleID % NelemX == 0)
                    {
                        qNodal[0] = 0.0;
                    }
                    if( globalEleID % NelemX == 1)
                    {
                        qNodal[0] = 1.0;
                    }


                }


                if (Testcase==33)
                {

                    for (int i=31; i<= 14881; i=i+150)
                    {
                        if( globalEleID == i)
                        {
                            qNodal[0] = 1.875;
                        }
                    }
                    for (int i=32; i<= 14882; i=i+150)
                    {
                        if( globalEleID == i)
                        {
                            qNodal[0] = 0.0;
                        }
                    }


                }

                if (Testcase==36)
                {

                    for (int i=31; i<= 14881; i=i+150)
                    {
                        if( globalEleID == i)
                        {
                            qNodal[0] = 1.875;
                        }
                    }
                    for (int i=32; i<= 14882; i=i+150)
                    {
                        if( globalEleID == i)
                        {
                            qNodal[0] = 0.1*PosPresTOL;
                        }
                    }


                }


                if (Testcase==43 ||Testcase==44 || Testcase==45 || Testcase == 46)
                {
                    if( globalEleID % 20 == 3)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID % 20 == 4)
                    {
                        qNodal[0] = 5.0;
                    }
                    if( globalEleID % 40 == 4)
                    {
                        qNodal[0] = 10.0;
                    }
                }
                if (Testcase == 47)
                {
                    if( globalEleID % 20 == 3)
                    {
                        qNodal[0] = 10.0;
                    }
                    if( globalEleID % 20 == 4)
                    {
                        qNodal[0] = 0.0;
                    }
                    if( globalEleID % 40 == 4)
                    {
                        qNodal[0] = 10.0;
                    }
                }

                q[qid] = qNodal[0];
                qid+=ngl2;   //first eq
                q[qid] = qNodal[1];
                qid+=ngl2;          //2nd eq
                q[qid] = qNodal[2];           //third eq
            }


        }
    }

}



void SW2D::InitQNodal(const dfloat x,const dfloat y, dfloat q[],const dfloat t,const dfloat b,const dfloat g_const, dfloat h_0)
{

    dfloat h,v,w;
    switch(Testcase)
    {

    case 1:      // periodic conv test
    {
        h=8.0+cos(x)*sin(y)*cos(t)-b;
        v= 0.5;
        w= 1.5;
        break;
    }

    case 2:      // WELL BALANCED (WITH BOTTOM)
    {

        h=10.0-b;


        v= 0.0;
        w= 0.0;
        break;
    }
    case 3:      // Entropy Glitch
    {
        if (x<0.0)
        {
            h=1.0-b;

        }
        else
        {
            h=0.1-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }
    case 4:      // WB _ smooth bottom
    {

        h=5.0-b;


        v= 0.0;
        w= 0.0;
        break;
    }

    case 5:      // Dry Lake + Lake at Rest
    {

        h=max(0.0,0.4-b);


        v= 0.0;
        w= 0.0;
        break;
    }

    case 6:      // Easy WB- debug test no bottom
    {

        h=5.0-b;


        v= 0.0;
        w= 0.0;
        break;
    }
    case 7:      // periodic conv test  _ NO BOTTOM
    {
        h=1.0+0.001*cos(2.0*PI*x)*sin(2.0*PI*y)*cos(t)-b;
        v= 0.5;
        w= 1.5;
        break;
    }
    case 8:      // periodic conv test  _ NO BOTTOM
    {
        h=20.0+cos(x)*sin(y)*cos(t)-b;
        v= 0.5;
        w= 1.5;
        break;
    }
    case 9:      // Dam Break (CARTESIAN)	(other direction)
    {
        if (y<0.0)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }

    case 20:      //  Dam Break (CARTESIAN)
    {
        if (x<0.0)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }
        v= 0.0;
        w= 0.0;
        break;
    }
    case 21:      // Dam Break (CARTESIAN)
    {
        if (x<0.0)
        {
            h=2.5-b;

        }
        else
        {
            h=1.5-b;
        }
        v= 0.0;
        w= 0.0;
        break;
    }
    case 29:      // Steeper Dam Break To Test Shock Capturing
    {
        if (x<0.0)
        {
            h=10.0-b;

        }
        else
        {
//        h=max(pow(10.0,-6),0.1-b);
            h=5.0;
        }


        v= 0.0;
        w= 0.0;
        break;
    }
    case 30:      // Steeper Dam Break To Test Shock Capturing
    {
        if (x<0.0)
        {
            h=10.0-b;

        }
        else
        {
//        h=max(pow(10.0,-6),0.1-b);
            h=0.0;
        }


        v= 0.0;
        w= 0.0;
        break;
    }
    case 28:      // Steeper Dam Break To Test Shock Capturing
    {
        if (x<0.5)
        {
            h=1.0-b;

        }
        else
        {
//        h=max(pow(10.0,-6),0.1-b);
            h=0.0;
        }


        v= 0.0;
        w= 0.0;
        break;
    }




    case 31:      // Two-dimensional oscillating lake  (Xing_PosPres paper, 6.8 or Gallardo for less mistakes)
    {
        //PARABOLIC BOWL

//    xL=-2.0
//    xR=2.0
//    yL=-2.0
//    yR=2.0
//    NelemX=100
//    NelemY=100
//    PeriodicBD_X=0
//    PeriodicBD_Y=0

        dfloat a=1.0;
        dfloat h0 = 0.1;
        dfloat sigma=0.5;
        dfloat omega = sqrt(2*g_const * h0 ) / a;
        h = max(0.0, sigma*h0 / pow(a,2) *(2 * x * cos(omega*t) + 2* y * sin(omega*t) - sigma) + h0 - b);

//    h = max(0.0, sigma*h0 / pow(a,2) *(2 * x * cos(omega*t) + 2* y * sin(omega*t) - sigma) + h0 - h0 * ( pow(x,2)+pow(y,2)/pow(a,2)));

//h=10.0-b;
        v= -sigma * omega * sin(omega*t);
        w= sigma * omega * cos(omega*t);


//    if (h>pow(10.0,-6)){
//        v= -sigma * omega * sin(omega*t);
//        w= sigma * omega * cos(omega*t);
//    }
//    else{
//        v =0.0;
//        w =0.0;
//    }

        break;
    }


    case 32:      // Flooding on a channel with three mounds
    {
        //(4.6) from Positivity-Preserving Well-Balanced Discontinuous Galerkin Methods for the Shallow Water Equations on unstructured triangular meshes
        h=0.0;
        v= 0.0;
        w= 0.0;


//    if (h>pow(10.0,-6)){
//        v= -sigma * omega * sin(omega*t);
//        w= sigma * omega * cos(omega*t);
//    }
//    else{
//        v =0.0;
//        w =0.0;
//    }

        break;
    }




    case 33:      // Dam Break over 3 Mounds
    {


//                    xL=0.0
//                    xR=75.0
//                    yL=0.0
//                    yR=30.0
//                    NelemX=150
//                    NelemY=100
//                    PeriodicBD_X=0
//                    PeriodicBD_Y=0


        if (x<16.0)
        {
            h=1.875-b;

        }
        else
        {
            h=0.0;
        }
        v= 0.0;
        w= 0.0;
        break;
    }



    case 34:      // 2D Solitary Wave Runup and Run-down
    {


//    xL=0.0
//    xR=25.0
//    yL=-15.0
//    yR=15.0

        dfloat A=0.064;
        dfloat x_c = 2.5;
        dfloat h0 = 0.32;
        dfloat gamma = sqrt((3.0*A)/(4*h0));
        dfloat displacement =  A/h0 /pow(cosh(gamma*(x-x_c)),2);
        dfloat zero=0.0;
        h= max (zero, h0 +   displacement  - b);
//	v= displacement * sqrt(g_const/h0);
        v= displacement * sqrt(g_const*h0);
        w= 0.0;
        break;
    }


    case 35:      // 1D Bowl
    {


//    xL=0.0
//    xR=25.0
//    yL=-15.0
//    yR=15.0
        dfloat B=5.0;
        dfloat Bterm = B*B/4.0/g_const;
        dfloat a=3000.0;
        dfloat h0 = 10.0;
        dfloat omega = sqrt(2*g_const * h0 ) / a;
        h = max(0.0,h0 - Bterm *g_const*(2*omega*t) - Bterm - B*x*0.5/a *sqrt(8*h0/g_const)*cos(omega*t)-b);

// end switch case


//    h = max(0.0, sigma*h0 / pow(a,2) *(2 * x * cos(omega*t) + 2* y * sin(omega*t) - sigma) + h0 - h0 * ( pow(x,2)+pow(y,2)/pow(a,2)));

//h=10.0-b;
        v= 0.0;
        w= 0.0;

        break;
    }


    case 36:      // Dam Break over 3 Mounds		// low minimal water height
    {


        if (x<16.0)
        {
            h=1.875-b;

        }
        else
        {
            h=0.1*PosPresTOL;
        }
        v= 0.0;
        w= 0.0;
        break;
    }


    case 43:      //curved dam break PARTIAL
    {
        if (x<1.0/25.0 * y*y - 0.25)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }

    case 44:      //curved dam break
    {
        if (x<1.0/25.0 * y*y - 0.25)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }


    case 45:      //curved dam NO BREAK
    {
        if (x<1.0/25.0 * y*y - 0.25)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }
    case 46:      //curved dam break PARTIAL
    {
        if (x<1.0/25.0 * y*y - 0.25)
        {
            h=10.0-b;

        }
        else
        {
            h=5.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }
    case 47:      //curved dam break PARTIAL, WET DRY
    {
        if (x<1.0/25.0 * y*y - 0.25)
        {
            h=10.0-b;

        }
        else
        {
            h=0.0-b;
        }

        v= 0.0;
        w= 0.0;
        break;
    }
    case 86:      // Completely wet WB Test for Ocean Mesh
    {
        h=20000.0+h_0 - b;
        v= 0.0;
        w= 0.0;
        break;
    }
    case 87:      // Wet/Dry Wellbalanced Test
    {
        dfloat zero = 0.0;
        h=max(zero,h_0-b);
        v= 0.0;
        w= 0.0;
        break;
    }
    case 88:      /// Conv Test for Ocean Mesh
    {
        //h=8.0+cos(2.0*PI*x)*sin(2.0*PI*y)*cos(t)-b;
        h= 20.0 + cos(0.001*x)*sin(0.001*y)*cos(t) - b;
        v= 0.5;
        w= 1.5;
        break;
    }

    case 89:      /// Conv Test for Ocean Mesh, for meters
    {
        //h=8.0+cos(2.0*PI*x)*sin(2.0*PI*y)*cos(t)-b;
        h= 20000.0 + cos(0.000001*x)*sin(0.000001*y)*cos(t) - b;
        v= 0.5;
        w= 1.5;
        break;
    }

    case 90:      /// Tsunami with okada  all wet (for meter mesh now)
    {
        //dfloat zero = 0.0;
        //h=max(zero,h_0-b);
        dfloat two = 2000.0;
        h=h_0+2000.0-b;
        v= 0.0;
        w= 0.0;
        break;
    }
    case 91:      /// Tsunami with okada  wet//dry
    {
        dfloat zero = 0.0;
        h=max(zero,h_0-b);
        v= 0.0;
        w= 0.0;
        break;
    }

    } //end switch


    q[0] = h;   //first eq
    q[1] = h*v;          //2nd eq
    q[2] = h*w;           //third eq

}

//    void InitB(const int Nelem,const int ngl,const int ngl2,const dfloat x[NoSpaceDofs],const dfloat y[NoSpaceDofs],dfloat b[NoSpaceDofs],const dfloat t){

void   SW2D::InitB(const int IsMeshSplit,const MeshPartitioning MeshSplit,const int Nelem,const int ngl,const int ngl2,const dfloat x[],const dfloat y[],dfloat * b)
{

    int globalEleID;
    for(int ie=0; ie<Nelem; ++ie)
    {
        if (IsMeshSplit==1)
        {
            globalEleID  = MeshSplit.MyElementLocalToGlobal[ie];
        }
        else
        {
            globalEleID  = ie+1;
        }
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                //erst alle Punkte die zu j=0 gehoeren, dann j=1, usw.
                int xid = ie*ngl2  + j*ngl+i;

                switch(Testcase)
                {
                case 1:      // convergence test
                {
                    b[xid] = 2.0+0.5*sin(2.0*PI*y[xid])+cos(2.0*PI*x[xid]);//+cos(x[xid])*sin(y[xid])*cos(t);
                    break;
                }

                case 2:      // WELL BALANCED (CARTESIAN 20x20)
                {
                    if ((globalEleID == 232)||(globalEleID==233) ||(globalEleID==231)||(globalEleID==212)||(globalEleID==213)||(globalEleID==211)||(globalEleID==252)||(globalEleID==251)||(globalEleID==253))
                    {
                        b[xid] =2.0+sin(y[xid])+cos(x[xid]);
                    }
                    else
                    {
                        b[xid] =0.0;//b[xid] =2.0+sin(y[xid])+cos(x[xid]);
                    }
                    break;
                }
                case 4:      // WELL BALANCED smooth bottom
                {
                    b[xid] =2.0+sin(y[xid])+cos(x[xid]);
                    break;
                }

                case 5:      // WELL BALANCED (CARTESIAN 4x4)
                {
                    b[xid] = 0.25 - 0.25*cos((2*x[xid]-1)*PI);
                    break;
                }
                case 8:      // convergence test
                {
                    b[xid] = 2.0+sin(y[xid])+cos(x[xid]);
                    break;
                }
                case 30:     // Steeper Dam Break To Test Shock Capturing
                {

//    if (x[xid] >= 2.5){
//        b[xid] = sin((x[xid]-2.5)/2.5 * PI);
//    }else{
//        b[xid] =0.0;
//    }
                    b[xid] =0.0;
                    break;
                }
                case 31:     // Two-dimensional oscillating lake  (Xing_PosPres paper, 6.8)
                {


                    dfloat a=1.0;
                    dfloat h0 = 0.1;
                    b[xid] = h0 * ( pow(x[xid],2)+pow(y[xid],2)/pow(a,2));

                    break;

                }
                case 32:     // Three Mound (4.6)
                {

                    dfloat Xm30Pow2 = pow((x[xid]-30),2);
                    dfloat Xm475Pow2 = pow((x[xid]-47.5),2);
                    dfloat Ym225Pow2 = pow((y[xid]-22.5),2);
                    dfloat Ym75Pow2 = pow((y[xid]-7.5),2);
                    dfloat Ym15Pow2 = pow((y[xid]-15),2);
                    dfloat m1 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym225Pow2  );
                    dfloat m2 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym75Pow2  );
                    dfloat m3 = 2.8-0.28 * sqrt(  Xm475Pow2 + Ym15Pow2  );
                    dfloat zero = 0.0;
                    b[xid] = max(max(zero,m1),max(m2,m3));


                    break;

                }

                case 33:     // Dam Break Three Mound (4.6)
                {

                    dfloat Xm30Pow2 = pow((x[xid]-30),2);
                    dfloat Xm475Pow2 = pow((x[xid]-47.5),2);
                    dfloat Ym225Pow2 = pow((y[xid]-22.5),2);
                    dfloat Ym75Pow2 = pow((y[xid]-7.5),2);
                    dfloat Ym15Pow2 = pow((y[xid]-15),2);
                    dfloat m1 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym225Pow2  );
                    dfloat m2 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym75Pow2  );
                    dfloat m3 = 2.8-0.28 * sqrt(  Xm475Pow2 + Ym15Pow2  );
                    dfloat zero = 0.0;
                    b[xid] = max(max(zero,m1),max(m2,m3));


                    break;

                }

                case 34:      // 2D Solitary Wave Runup and Run-down
                {
                    dfloat r_c=3.6;
                    dfloat x_c = 12.5;
                    dfloat y_c = 15;
                    dfloat r = sqrt(pow(x[xid]-x_c,2) + pow(y[xid]-y_c,2));
                    if (r <= r_c)
                    {
                        b[xid] = 0.93*(1.0-r/r_c);
                    }
                    else
                    {
                        b[xid] = 0.0;
                    }

                    break;
                }
                case 35:      // 1D Bowl
                {


//    xL=0.0
//    xR=25.0
//    yL=-15.0
//    yR=15.0
                    b[xid] = 10.0*x[xid]*x[xid]/9000000.0;



                    break;
                }

                case 36:     // Dam Break Three Mound (4.6)
                {

                    dfloat Xm30Pow2 = pow((x[xid]-30),2);
                    dfloat Xm475Pow2 = pow((x[xid]-47.5),2);
                    dfloat Ym225Pow2 = pow((y[xid]-22.5),2);
                    dfloat Ym75Pow2 = pow((y[xid]-7.5),2);
                    dfloat Ym15Pow2 = pow((y[xid]-15),2);
                    dfloat m1 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym225Pow2  );
                    dfloat m2 = 1.0-0.1 * sqrt(  Xm30Pow2 + Ym75Pow2  );
                    dfloat m3 = 2.8-0.28 * sqrt(  Xm475Pow2 + Ym15Pow2  );
                    dfloat zero = 0.0;
                    b[xid] = max(max(zero,m1),max(m2,m3));


                    break;

                }

                case 43:      //curved dam break PARTIAL
                {
                    if (x[xid] >= 2.25)
                    {
                        b[xid] = 2.0 + log(x[xid]-1.25);

                    }
                    else
                    {
                        b[xid] =0.0;
                    }
                    if ((x[xid]==2.25)&&( globalEleID%10 ==9))
                    {
                        b[xid] =0.0;
                    }
                    break;
                }

                case 44:     //curved dam break
                {
                    if (x[xid] >= 2.25)
                    {
                        b[xid] = 2.0 + log(x[xid]-1.25);

                    }
                    else
                    {
                        b[xid] =0.0;
                    }
                    if ((x[xid]==2.25)&&( globalEleID%10 ==9))
                    {
                        b[xid] =0.0;
                    }
                    break;
                }


                case 45:      //curved dam NO BREAK
                {
                    if (x[xid] >= 2.25)
                    {
                        b[xid] = 2.0 + log(x[xid]-1.25);

                    }
                    else
                    {
                        b[xid] =0.0;
                    }
                    if ((x[xid]==2.25)&&( globalEleID%10 ==9))
                    {
                        b[xid] =0.0;
                    }
                    break;
                }
                default:
                {
                    b[xid] =0.0;
                    break;
                }
                }



//            b[xid] =2.0+sin(y[xid])+cos(x[xid]);//+cos(x[xid])*sin(y[xid])*cos(t);
//            b[xid]=0.0;

            }


        }
    }

}

void  SW2D::CalcBDerivatives(const int Nelem,const int ngl,const int ngl2,const dfloat g_const,const dfloat x[],const dfloat y[],const dfloat * b,const dfloat D[],const dfloat y_eta[],const dfloat y_xi[],const dfloat x_eta[],const dfloat x_xi[],dfloat Bx[],dfloat By[],const dfloat J[])
{


//CALC BOTTOM DERIVATIVES
    for (int ie=0; ie<Nelem; ie++)
    {
        for (int j=0; j<ngl; ++j)
        {
            for (int i=0; i<ngl; ++i)
            {

                int id = ie*ngl2  ; //ankerpunkt

                //bottom terms for EQ2
                dfloat B_xi=0.0;
                dfloat BYeta_xi=0.0;
                dfloat BXxi_eta=0.0;
                //bottom terms for EQ3
                dfloat B_eta=0.0;
                dfloat BYxi_eta=0.0;
                dfloat BXeta_xi=0.0;


                for (int l=0; l<ngl; ++l)
                {
                    int il=i*ngl+l;
                    int jl=j*ngl+l;
                    int li=l*ngl+i;
                    int lj=l*ngl+j;


                    // B_x
                    B_xi+=D[il] * b[id+jl];
                    BYeta_xi+=D[il] * b[id+jl]* y_eta[id+jl];
                    BXeta_xi+=D[il] * b[id+jl]* x_eta[id+jl];




                    //B_y
                    B_eta+=D[jl] * b[id+ li];
                    BYxi_eta+=D[jl] * b[id+ li]*y_xi[id+li];
                    BXxi_eta+=D[jl] * b[id+ li]*x_xi[id+li];
                }


                int xid=id +j*ngl+i;

                Bx[xid]=  0.5*J[xid]*(y_eta[xid]*B_xi + BYeta_xi - y_xi[xid]* B_eta - BYxi_eta);
                By[xid]=  0.5*J[xid]*(-x_eta[xid]*B_xi - BXeta_xi + x_xi[xid]* B_eta + BXxi_eta);


            }
        }
    }

}



