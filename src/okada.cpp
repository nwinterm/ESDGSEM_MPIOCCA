
#include "okada.h"
//#include "Constants.hpp"





okada::okada(const int Nelem)
{
NumElements = Nelem;

}




void okada::okadamapFull(const int ngl, const dfloat x[],const dfloat y[], dfloat q[])
{
    /*!
     * return displacement dZ for a surface displacement at (xloc, yloc)
     * given okadaparams
     * adopted from okada.py
     */
    const dfloat zero	= 0.0;
    const dfloat osixty= 1.0/60.0; //0.016666666667;
    const dfloat rad	= M_PI/180. ;//0.01745329252;
    const dfloat earth_radius= 6.378e6;



    for(int fault = 1; fault <=1 ; ++fault)
    {

        // WHAT DO THESE DO?!



        char dummy[BUFSIZ];
        sprintf(dummy, "%s%d.cfg", "./data/indianOceanUSGS", fault);
        dfloat w, l, rd, dl, th, d, y0, x0, hh;

        std::ifstream usgsfile(dummy);
        usgsfile >> x0; // epicenter longitude
        usgsfile >> y0; // epicenter latitude
        usgsfile >> hh; // depth
        usgsfile >> l; // fault_length
        usgsfile >> w; // fault_width

        usgsfile >> th; // strike angle
        usgsfile >> rd; // slip angle
        usgsfile >> dl; // dip angle

        usgsfile >> d; // dislocation

        usgsfile.close();

        /// get epicenter coordinates in meters
        const dfloat x0m = (x0/earth_radius)*(180./M_PI);
        const dfloat y0m = (atan(sinh(y0/earth_radius))) * (180./M_PI);

        const dfloat xlower = x0m-l;
        const dfloat xupper = x0m+l;
        const dfloat ylower = y0m-w;
        const dfloat yupper = y0m+w;

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

        //const dfloat xl = rr*cos(rad*yo)*(x0-xo)*rad + del_x;   ///original source from geoclaw, okada.py
        const dfloat yl = earth_radius*(y0-yo)*rad - del_y; ///from gandham, dont know why

        for(int ie=0; ie<NumElements; ++ie) /// loop through all elements
        {
            for(int i=0; i<ngl; ++i)  /// loop through all nodes of element
            {
                for(int j=0; j<ngl; ++j)  /// loop through all nodes of element
                {
                    const int ngl2=ngl*ngl;
                    const int Neq=3;
                    int qid = ie*ngl2*Neq   +j*ngl+i;
                    int xid = ie*ngl2  + j*ngl+i;
                    // get actual x/y in meters
                    dfloat xloc = x[xid];
                    dfloat yloc = y[xid];

                    // check if this is in range of okada model
                    if(xloc >= xlower && xloc <= xupper &&
                            yloc >= ylower && yloc <= yupper )
                    {
                        /// convert to degrees
                        xloc = (xloc/earth_radius)*(180./M_PI);
                        yloc = (atan(sinh(yloc/earth_radius))) * (180./M_PI);

                        const dfloat xl = earth_radius*cos(rad*yloc)*(x0-xo)*rad + del_x;
                        const dfloat yy = earth_radius*(yloc-yo)*rad;
                        const dfloat xx = earth_radius*cos(rad*yloc)*(xloc-xo)*rad;


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
                        q[qid] += (us+ud);
                    }
                }
            }
        }
    }


}


dfloat strike_slip (const dfloat x1,
                    const dfloat x2,
                    const dfloat x3,
                    const dfloat y1,
                    const dfloat y2,
                    const dfloat dp,
                    const dfloat dd)
{
    /*!
     * Used for Okada's model
     * Code borrowed from gandham who borrowed from okada.py in geoclaw
     */
    const dfloat sn = sin(dp);
    const dfloat cs = cos(dp);
    const dfloat p = x2*cs + dd*sn;
    const dfloat q = x2*sn - dd*cs;
    const dfloat d_bar = y2*sn - q*cs;
    const dfloat r = sqrt(y1*y1 + y2*y2 + q*q);
    const dfloat xx = sqrt(y1*y1 + q*q);
    const dfloat a4 = 0.5*1./cs*(log(r+d_bar) - sn*log(r+y2));
    const dfloat f =
        -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*M_PI);

    return f;
}


dfloat dip_slip (const dfloat x1,
                 const dfloat x2,
                 const dfloat x3,
                 const dfloat y1,
                 const dfloat y2,
                 const dfloat dp,
                 const dfloat dd)
{
    /*!
     * Based on Okada's paper (1985)
     * Code borrowed from gandham who  borrowed from okada.py in geoclaw
     */
    const dfloat sn = sin(dp);
    const dfloat cs = cos(dp);

    const dfloat p = x2*cs + dd*sn;
    const dfloat q = x2*sn - dd*cs;
    const dfloat d_bar = y2*sn - q*cs;
    const dfloat r = sqrt(y1*y1 + y2*y2 + q*q);
    const dfloat xx = sqrt(y1*y1 + q*q);
    const dfloat a5 = 0.5*2/cs*atan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs);
    const dfloat f = -(d_bar*q/r/(r+y1) + sn*atan(y1*y2/q/r) - a5*sn*cs)/(2.0*M_PI);

    return f;
}
