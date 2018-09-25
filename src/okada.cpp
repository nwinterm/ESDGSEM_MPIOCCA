
#include "okada.h"
//#include "Constants.hpp"





okada::okada(const int tcase, const dfloat postol)
{
    Testcase = tcase;
    PosPresTOL = postol;

}


void okada :: InitQ(const int IsMeshSplit,const MeshPartitioning MeshSplit,const int Nelem,const int ngl,const int ngl2,
                    const dfloat x[],const dfloat y[], dfloat q[],const dfloat t,
                    const dfloat * b, const dfloat g_const, dfloat h_0)
{


}


/*!
 * borrowed code from Gandham who borrowed from Frank X Giraldo's DGCOM
 */
void okada::set_ics_from_okada_model()
{



    datafloat xlower, xupper, ylower, yupper;
    int xsize, ysize;
    fmatrix dZ, XGRID, YGRID;

    okadamapFull( xsize, ysize, xlower,
                 xupper, ylower, yupper);

    return;


}







void okada::okadamapFull(int &xsize, int &ysize,
                            dfloat &xlower,
                            dfloat &xupper,
                            dfloat &ylower,
                            dfloat &yupper)
{
    /*!
     * return displacement dZ for a surface displacement at (xloc, yloc)
     * given okadaparams
     * adopted from okada.py
     */
    const dfloat zero	= 0.0;
    const dfloat osixty= 1.0/60.0; //0.016666666667;
    const dfloat rad	= M_PI/180. ;//0.01745329252;
    const dfloat rr	= 6.378e6;



    for(int fault = 1; fault <=5 ; ++fault)
    {

        char dummy[BUFSIZ];
        sprintf(dummy, "%s%d.cfg", "./data/indianOceanUSGS", fault);
        dfloat w, l, rd, dl, th, d, y0, x0, hh;

        std::ifstream usgsfile(dummy);

        usgsfile >> dummy;
        usgsfile >> w; // fault_width
        usgsfile >> dummy;
        usgsfile >> l; // fault_length
        usgsfile >> dummy;
        usgsfile >> rd; // slip angle
        usgsfile >> dummy;
        usgsfile >> dl; // dip angle
        usgsfile >> dummy;
        usgsfile >> th; // strike direction
        usgsfile >> dummy;
        usgsfile >> d; // dislocation
        usgsfile >> dummy;
        usgsfile >> y0; // epicenter latitude
        usgsfile >> dummy;
        usgsfile >> x0; // epicenter longitude
        usgsfile >> dummy;
        usgsfile >> hh; // focal depth
        usgsfile >> dummy;
        usgsfile >> xsize; // grid points in x-dir
        usgsfile >> dummy;
        usgsfile >> ysize; // grid points in y-dir
        usgsfile >> dummy;
        usgsfile >> ylower;
        usgsfile >> dummy;
        usgsfile >> yupper;
        usgsfile >> dummy;
        usgsfile >> xlower;
        usgsfile >> dummy;
        usgsfile >> xupper;

        usgsfile.close();

        const dfloat xo = xlower;//XGRID(1,1);//xlower;
        const dfloat yo = yupper;//YGRID(1,1);//yupper;

        const dfloat ang_l = rad*dl;
        const dfloat ang_r = rad*rd;
        const dfloat ang_t = rad*th;
        const dfloat halfl = 0.5*l;

        /// Calculate focal depth used for Okada's model
        hh = hh+w*sin(ang_l);

        ///   displacement due to different epicenter definition
        const dfloat del_x = w*cos(ang_l)*cos(ang_t);
        const dfloat del_y = w*cos(ang_l)*sin(ang_t);

        const dfloat ds = d*cos(ang_r);
        const dfloat dd = d*sin(ang_r);
        const dfloat sn = sin(ang_l);
        const dfloat cs = cos(ang_l);

        //  const dfloat xl = rr*cos(rad*yo)*(x0-xo)*rad + del_x;
        const dfloat yl = rr*(y0-yo)*rad - del_y;

        int countInterior = 0;
        for(int k=1; k<=K; ++k)
        {
            for(int n=1; n<=Np; ++n)
            {
                // get actual x/y in meters
                dfloat xloc = x(n,k);
                dfloat yloc = y(n,k);

                /// convert to degrees
                xloc = (xloc/earth_radius)*(180./M_PI);
                yloc = (atan(sinh(yloc/earth_radius))) * (180./M_PI);

                // check if this is in range of okada model
                if(xloc >= xlower && xloc <= xupper &&
                        yloc >= ylower && yloc <= yupper )
                {

                    countInterior++;
                    const dfloat xl = rr*cos(rad*yloc)*(x0-xo)*rad + del_x;
                    const dfloat yy = rr*(yloc-yo)*rad;
                    const dfloat xx = rr*cos(rad*yloc)*(xloc-xo)*rad;


                    dfloat x1 = (xx-xl)*sin(ang_t)+(yy-yl)*cos(ang_t);
                    dfloat x2 = (xx-xl)*cos(ang_t)-(yy-yl)*sin(ang_t);
                    x2 = -x2;
                    dfloat x3 = zero;

                    const dfloat p  = x2*cs+hh*sn;

                    dfloat f1 = strike_slip (x1,x2,x3,x1+halfl,p,ang_l,hh);
                    dfloat f2 = strike_slip (x1,x2,x3,x1+halfl,p-w,ang_l,hh);
                    dfloat f3 = strike_slip (x1,x2,x3,x1-halfl,p,ang_l,hh);
                    dfloat f4 = strike_slip (x1,x2,x3,x1-halfl,p-w,ang_l,hh);

                    dfloat g1 = dip_slip (x1,x2,x3,x1+halfl,p,ang_l,hh);
                    dfloat g2 = dip_slip (x1,x2,x3,x1+halfl,p-w,ang_l,hh);
                    dfloat g3 = dip_slip (x1,x2,x3,x1-halfl,p,ang_l,hh);
                    dfloat g4 = dip_slip (x1,x2,x3,x1-halfl,p-w,ang_l,hh);


                    const dfloat us = (f1-f2-f3+f4)*ds;
                    const dfloat ud = (g1-g2-g3+g4)*dd;

                    // This seems to update the current water height by strike/dip slips ?? displacements?!
                    Q(n,k) += (us+ud);
                }
            }
        }
    }

    dfloat hmax = -1e9;
    dfloat hmin = 1e9;

    for(int k=1; k<=K; ++k)
    {
        for(int n=1; n<=Np; ++n)
        {

            dfloat pe = Q(n,k);
            dfloat pp = pe - B(n,k);


            if(pp < h_thresh)
                pp = h_thresh;


            Q(n,k) = pp + B(n,k);
            Q(n+Npad, k) = 0.;
            Q(n+2*Npad, k) = 0.;

            //      assert(Q(n,k) > 0);

            if(hmax < Q(n,k))
                hmax = Q(n,k);

            if(hmin > Q(n,k))
                hmin = Q(n,k);
        }

        for(int fld=1; fld<=Nfields; ++fld)
        {
            for(int n=1; n<=Ngauss*Nfaces; ++n)
            {

                dfloat gQn = 0.;

                for(int m=1; m<=Np; ++m)
                {

                    gQn += gInterp(n,m)*Q(m + (fld-1)*Npad, k);
                }

                gQ(n + (fld-1)*Ngauss_NfacesPad, k) = gQn;

                assert(gQn == gQn);

            }
        }
    }

    std::cout << " hmax = " << hmax << std::endl;
    std::cout << " hmin = " << hmin << std::endl;
}
