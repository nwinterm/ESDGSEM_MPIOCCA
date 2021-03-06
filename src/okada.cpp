//==================================================================================================================================
// Copyright (c) 2019 Niklas Wintermeyer
// Copyright (c) 2019 Gregor Gassner
// Copyright (c) 2019 Andrew Winters
//
// This file is part of ESDGSEM_MPIOCCA (github.com/ESDGSEM_MPIOCCA). ESDGSEM_MPIOCCA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version.
//
// ESDGSEM_MPIOCCA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
//
// You should have received a copy of the GNU General Public License along with ESDGSEM_MPIOCCA. If not, see <http://www.gnu.org/licenses/>.



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
    const dfloat earth_radius= 6.378e6;	// was previously = 6.371e6 but that is not correct for lon/lat measurements



    for(int fault = 1; fault <=5 ; ++fault)
    {

        // WHAT DO THESE DO?!



        char dummy[BUFSIZ];
        sprintf(dummy, "%s%d.txt", "./data/indianOceanFault", fault);
        dfloat w, l, rd, dl, th, d, y0, x0, hh;

        std::ifstream usgsfile(dummy);

        if (!usgsfile)
        {
            std::string error_message("ERROR: Okada Input file not found: ");
            error_message += dummy;
            throw std::invalid_argument(error_message);
        }
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

//        cout << "Epicenter of okada earthquake: " << x0 << ", " << y0 << "\n";
//        cout << "depth of okada earthquake: " << hh << "\n";
//        cout << "fault_length and fault_width of okada earthquake: " << l << ", " << w << "\n";
//        cout << "Strike Slip and Dip angles of okada earthquake: " << th << ", " << rd << ", " << dl << "\n";
//        cout << "dislocation of okada earthquake: " << d << "\n";

        /// get epicenter coordinates in meters
        const dfloat x0m = x0*earth_radius*M_PI/180.0;
        const dfloat y0m = earth_radius*asinh(tan(y0/180.0*M_PI));

//        cout << "Epicenter (KM) of okada earthquake: " << x0m << ", " << y0m << "\n";


        const dfloat xlower = x0m-l;
        const dfloat xupper = x0m+l;
        const dfloat ylower = y0m-w;
        const dfloat yupper = y0m+w;

        // get lower values in degrees
        const dfloat xo = (xlower/earth_radius)*(180./M_PI);
        const dfloat yo = (atan(sinh(ylower/earth_radius))) * (180./M_PI);




        const dfloat ang_l = rad*dl;
        const dfloat ang_r = rad*rd;
        const dfloat ang_t = rad*th;
        const dfloat halfl = 0.5*l;

        /// Calculate focal depth used for Okada's model
        hh = hh+w*sin(ang_l);

        ///   displacement due to different epicenter definition
        const dfloat del_x = w*cos(ang_l)*cos(ang_t);
        const dfloat del_y = w*cos(ang_l)*sin(ang_t);


        dfloat xl = earth_radius*cos(rad*yo)*(x0-xo)*rad + del_x;   		///original source from geoclaw, okada.py
        const dfloat yl = earth_radius*(y0-yo)*rad - del_y; 			///from gandham, dont know why
        const dfloat ds = d*cos(ang_r);
        const dfloat dd = d*sin(ang_r);
        const dfloat sn = sin(ang_l);
        const dfloat cs = cos(ang_l);





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
                    // get actual x/y in METERS
                    dfloat xloc = x[xid]*1000.0;
                    dfloat yloc = y[xid]*1000.0;





                    // check if this is in range of okada model
                    if(xloc >= xlower && xloc <= xupper &&
                            yloc >= ylower && yloc <= yupper )
                    {
                        /// convert to degrees
                        xloc = (xloc/earth_radius)*(180./M_PI);
                        yloc = (atan(sinh(yloc/earth_radius))) * (180./M_PI);


                        ///#!-----added by Xiaoming Wang, Nov 11 2006----
                        xl = earth_radius*cos(rad*yloc)*(x0-xo)*rad + del_x ; ///# TAL - Fixed sign, 9/07
                        ///#!---------------------------------------------

                        // const dfloat xl = earth_radius*cos(rad*yloc)*(x0-xo)*rad + del_x;
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
                        q[qid] += (us+ud)/1000.0;
                    }
                }
            }
        }
    }


}


dfloat okada:: strike_slip (const dfloat x1,
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


dfloat okada:: dip_slip (const dfloat x1,
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
