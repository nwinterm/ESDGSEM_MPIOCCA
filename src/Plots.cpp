#include "plots.h"



void PlotSolution(const int Nelem, const int ngl,const int PlotVar, const dfloat x[], const dfloat y[], const dfloat Q[], const dfloat b[], const int plotCount, const dfloat h_0)
{

    if (PlotVar>0)
    {
        ostringstream os;
        os << "movie/plot_" << plotCount << ".tec";
        string fName = os.str();


        ofstream plotfile;
        plotfile.open (fName.c_str());


        int ngl2=ngl*ngl;
        dfloat Qinv;
        switch(PlotVar)
        {
        case 1:
        {

            plotfile <<"TITLE = H_solution.tec\n";
            plotfile <<"VARIABLES = \"x\",\"y\",\"H\",\"hu\",\"hv\",\"bottom\"\n";
            for (int ie=0; ie<Nelem; ie++)
            {
                plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
                for(int j=0; j<ngl; ++j)
                {
                    for(int i=0; i<ngl; ++i)
                    {
                        int id = ie*ngl2*Neq   +j*ngl+i;
                        int xid = ie*ngl2   +j*ngl+i;
                        dfloat H;
                        if (Q[id] > 0.0)
                        {
                            Qinv = 1.0/Q[id];
                        }
                        else
                        {
                            Qinv = 0.0;
                        }


                        H=Q[id]+b[xid];

                        plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]<<" " << Q[id+ngl2+ngl2]<<" "<<b[xid]<<" \n";
                        //plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]*Qinv<<" " << Q[id+ngl2+ngl2]*Qinv<<" "<<b[xid]<<" \n";

                    }
                }
            }
            break;
        }
        case 2:
        {

            plotfile <<"TITLE = H_solution.tec\n";
            plotfile <<"VARIABLES = \"x\",\"y\",\"FreeSurface\",\"hu\",\"hv\",\"bottom\"\n";
            for (int ie=0; ie<Nelem; ie++)
            {
                plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
                for(int j=0; j<ngl; ++j)
                {
                    for(int i=0; i<ngl; ++i)
                    {
                        int id = ie*ngl2*Neq   +j*ngl+i;
                        int xid = ie*ngl2   +j*ngl+i;
                        dfloat H;


                        H=Q[id]+b[xid]-h_0;

                        plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]<<" " << Q[id+ngl2+ngl2]<<" "<<b[xid]<<" \n";
                        //plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]*Qinv<<" " << Q[id+ngl2+ngl2]*Qinv<<" "<<b[xid]<<" \n";

                    }
                }
            }
            break;
        }
        }

        plotfile.close();

    }
}

void PlotFriction(const int Nelem, const int ngl,const int PlotVar, const dfloat x[], const dfloat y[], const dfloat Friction[], const int plotCount)
{

    ostringstream os;
    os << "movie/Friction_" << plotCount << ".tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());


    int ngl2=ngl*ngl;

    plotfile <<"TITLE = Friction.tec\n";
    plotfile <<"VARIABLES = \"x\",\"y\",\"Friction1\",\"Friction2\"\n";
    for (int ie=0; ie<Nelem; ie++)
    {
        plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*(Neq-1)   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;



                plotfile <<x[xid]<<" "<<y[xid]<<" " << Friction[id]<<" " << Friction[id+ngl2]<<" \n";
                //plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]*Qinv<<" " << Q[id+ngl2+ngl2]*Qinv<<" "<<b[xid]<<" \n";

            }
        }
    }


    plotfile.close();


}


void PlotArrivalTimings(const int Nelem, const int ngl, const dfloat x[], const dfloat y[], const dfloat Arrivaltimes[])
{

    ostringstream os;
    os << "movie/Arrivals.tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());


    int ngl2=ngl*ngl;

    plotfile <<"TITLE = Arrivals.tec\n";
    plotfile <<"VARIABLES = \"x\",\"y\",\"Arrivaltime\",\n";
    for (int ie=0; ie<Nelem; ie++)
    {
        plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int xid = ie*ngl2   +j*ngl+i;



                plotfile <<x[xid]<<" "<<y[xid]<<" " << Arrivaltimes[xid] <<" \n";
                //plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]*Qinv<<" " << Q[id+ngl2+ngl2]*Qinv<<" "<<b[xid]<<" \n";

            }
        }
    }


    plotfile.close();


}

void PlotSolutionWithExact(const int Nelem, const int ngl,const int PlotVar, const dfloat x[], const dfloat y[], const dfloat Q[], const dfloat b[], const int plotCount, const dfloat qExact[])
{

    ostringstream os;
    os << "movie/plot_" << plotCount << ".tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());


    int ngl2=ngl*ngl;
    dfloat Qinv;
    switch(PlotVar)
    {
    case 1:

        plotfile <<"TITLE = H_solution.tec\n";
        plotfile <<"VARIABLES = \"x\",\"y\",\"H\",\"u\",\"v\",\"bottom\",\"exact\"\n";

        for (int ie=0; ie<Nelem; ie++)
        {
            plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
            for(int j=0; j<ngl; ++j)
            {
                for(int i=0; i<ngl; ++i)
                {
                    int id = ie*ngl2*Neq   +j*ngl+i;
                    int xid = ie*ngl2   +j*ngl+i;
                    dfloat H;
                    if (Q[id] > 0.0)
                    {
                        Qinv = 1.0/Q[id];
                    }
                    else
                    {
                        Qinv = 0.0;
                    }

                    if (Q[id] > pow(10.0,2))
                    {
                        H=-1;
                    }
                    else
                    {
                        H=Q[id]+b[xid];
                    }

                    plotfile <<x[xid]<<" "<<y[xid]<<" "<<H<< " " << Q[id+ngl2]*Qinv<<" " << Q[id+ngl2+ngl2]*Qinv<<" "<<b[xid]<< " " << qExact[id]+b[xid] << " \n";



                    //plotfile <<x[xid]<<" "<<y[xid]<<" "<<Q[id]+b[xid]<< " " << Q[id+ngl2]<<" " << Q[id+ngl2+ngl2]<<" "<<b[xid]<<" \n";
                }
            }
        }
        break;
    }

    plotfile.close();


}

void PlotViscosity(const int Nelem, const int ngl,const int PlotVar, const dfloat x[], const dfloat y[], const dfloat Qx[], const dfloat Qy[], const int plotCount)
{

    ostringstream os;
    os << "movie/gradients_" << plotCount << ".tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());


    int ngl2=ngl*ngl;




    plotfile <<"TITLE = H_solution.tec\n";
    plotfile <<"VARIABLES = \"x\",\"y\",\"Qx1\",\"Qx2\",\"Qx3\",\"Qy1\",\"Qy2\",\"Qy3\"\n";

    for (int ie=0; ie<Nelem; ie++)
    {
        plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;
                plotfile <<x[xid]<<" "<<y[xid]<<" "<<Qx[id]<<" "<<Qx[id+ngl2]<<" "<<Qx[id+ngl2+ngl2]<<" "<<Qy[id]<<" "<<Qy[id+ngl2]<<" "<<Qy[id+ngl2+ngl2]<<" \n";


            }
        }
    }



    plotfile << "Test";
    plotfile.close();


}


void PlotViscoseParameter(const int Nelem, const int ngl, const dfloat x[], const dfloat y[], const dfloat Eps[], const int plotCount)
{

    ostringstream os;
    os << "movie/visc_" << plotCount << ".tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());


    int ngl2=ngl*ngl;


    plotfile <<"TITLE = ArtVisc.tec\n";
    plotfile <<"VARIABLES = \"x\",\"y\",\"Eps\"\n";

    for (int ie=0; ie<Nelem; ie++)
    {
        plotfile <<"ZONE I ="<<ngl<<",J="<<ngl<<",F=POINT\n";
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int id = ie*ngl2*Neq   +j*ngl+i;
                int xid = ie*ngl2   +j*ngl+i;
                plotfile <<x[xid]<<" "<<y[xid]<<" "<<Eps[ie]<<" \n";


            }
        }
    }









    plotfile << "Test";
    plotfile.close();





}


void PlotEntropy(const int NumPlots, const dfloat EntropyTimes[], const dfloat TotalEntropy[])
{

    ostringstream os;
    os << "movie/Entropyplot.tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());




    plotfile <<"TITLE = Entropyplot.tec\n";
    plotfile <<"VARIABLES = \"t\",\"Entropy\"\n";


    for(int i=0; i<NumPlots; ++i)
    {
        plotfile <<EntropyTimes[i]<<" "<<TotalEntropy[i]<<" \n";
    }







    plotfile.close();





}

void PlotMass(const int NumPlots, const dfloat PlotTimes[], const dfloat TotalMass[])
{

    ostringstream os;
    os << "movie/Massplot.tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());




    plotfile <<"TITLE = Massplot.tec\n";
    plotfile <<"VARIABLES = \"t\",\"Massdifference\"\n";


    for(int i=0; i<NumPlots; ++i)
    {
        plotfile <<PlotTimes[i]<<" "<<(TotalMass[i]-TotalMass[0])/TotalMass[0]<<" \n";
    }







    plotfile.close();





}




void PlotTimeSeries(const int NumPlots, const dfloat TimeSeriesTimes[],
                     const dfloat ChennaiTimeSeries[],const dfloat TuticorinTimeSeries[],
                     const dfloat VisakhapatnamTimeSeries[],const dfloat ParadipTimeSeries[],
                     const dfloat KochiTimeSeries,const dfloat MormugaoTimeSeries,
                     const dfloat OkhaTimeSeries)
{

    ostringstream os;
    os << "movie/Timeseries.tec";
    string fName = os.str();


    ofstream plotfile;
    plotfile.open (fName.c_str());




    plotfile <<"TITLE = Timeseries.tec\n";
    plotfile <<"VARIABLES = \"t\",\"Chennai\",\"Tuticorin\",\"Visakhapatnam\",\"Paradip\",\"Kochi\",\"Mormugao\",\"Okha\"\n";


    for(int i=0; i<NumPlots; ++i)
    {
        plotfile <<TimeSeriesTimes[i]<<" "<<ChennaiTimeSeries[i]<<" "<<TuticorinTimeSeries[i]<<" "<<VisakhapatnamTimeSeries[i]<<" "<<ParadipTimeSeries[i]<<" "<<KochiTimeSeries[i]<<" "<<MormugaoTimeSeries[i]<<" "<<OkhaTimeSeries[i]<<" \n";
    }







    plotfile.close();





}
