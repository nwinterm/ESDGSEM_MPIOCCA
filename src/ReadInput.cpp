#include "ReadInput.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <string>
using namespace std;


string getStrValue(string line)
{
    string returnValue;
    int pos=line.find("=");

//    std::string value = current_string.substr(pos+1);


    std::stringstream ssValue(line.substr(pos+1));
    ssValue>>returnValue;

    return returnValue;


}

double getDoubleValue(string line)
{
    double returnValue;
    int pos=line.find("=");

//    std::string value = current_string.substr(pos+1);


    std::stringstream ssValue(line.substr(pos+1));
    ssValue>>returnValue;

    return returnValue;
}

dfloat getDfloatValue(string line)
{
    dfloat returnValue;
    int pos=line.find("=");

//    std::string value = current_string.substr(pos+1);


    std::stringstream ssValue(line.substr(pos+1));
    ssValue>>returnValue;

    return returnValue;
}


int getIntValue(string line)
{
    int returnValue;
    int pos=line.find("=");

//    std::string value = current_string.substr(pos+1);


    std::stringstream ssValue(line.substr(pos+1));
    ssValue>>returnValue;

    return returnValue;
}

bool getBoolValue(string line)
{
    bool returnValue;
    int pos=line.find("=");

//    std::string value = current_string.substr(pos+1);


    std::stringstream ssValue(line.substr(pos+1));
    ssValue>>returnValue;

    return returnValue;
}


void ReadInputFile(int *N,
                   string *meshFile,
                   dfloat *CFL,
                   dfloat *DFL,
                   dfloat *T,
                   dfloat *g_const,
                   int *ArtificialViscosity,
                   int *PositivityPreserving,
                   dfloat *PosPresTOL,
                   dfloat *epsilon_0,
                   dfloat *sigma_min,
                   dfloat *sigma_max,
                   int *PlotVar,
                   int *EntropyPlot,
                   int *NumPlots,
                   int *NumTimeChecks,
                   int *Testcase,
                   int *ES,
                   int *NumFlux,
                   int *FluxDifferencing,
                   int *Cartesian,
                   int *rkorder,
                   int *rkSSP,
                   int *NEpad,
                   int *NEsurfpad,
                   int *Nedgepad,
                   int *NAvgPad,
                   int *KernelVersion,
                   int * KernelVersionSTD,
                   int * DiscBottom,
                   int * ReadInBottom,
                   int * PartialDryTreatment,
                   int * FrictionTerms,
                   int * ConvertToKM)
{

//,dfloat T, dfloat g_const


    std::ifstream InputStream;
    string filename="input.dat";
    InputStream.open(filename.c_str());

    if (!InputStream)
    {
        std::string error_message("ERROR: Input file not found: ");
        error_message += filename;
        throw std::invalid_argument(error_message);
    }


    std::string current_string;

    std::getline(InputStream, current_string);
    *meshFile= getStrValue(current_string);
    std::getline(InputStream, current_string);
    if (*N ==0)
    {
        *N= getIntValue(current_string);
    }

    std::getline(InputStream, current_string);
    *CFL= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *DFL= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *T= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *g_const= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *ArtificialViscosity= getIntValue(current_string);

    std::getline(InputStream, current_string);
    *PositivityPreserving= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *PosPresTOL= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *epsilon_0= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *sigma_min= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *sigma_max= getDfloatValue(current_string);
    std::getline(InputStream, current_string);
    *PlotVar= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *EntropyPlot= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *NumPlots= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *NumTimeChecks= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *Testcase= getIntValue(current_string);

    std::getline(InputStream, current_string);
    *ES= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *NumFlux= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *FluxDifferencing= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *Cartesian= getBoolValue(current_string);

    std::getline(InputStream, current_string);
    *rkorder= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *rkSSP= getIntValue(current_string);
    std::getline(InputStream, current_string);
    if (*NEpad==0)
    {
        *NEpad= getIntValue(current_string);
    }

    std::getline(InputStream, current_string);
    *NEsurfpad= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *Nedgepad= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *NAvgPad= getIntValue(current_string);
    std::getline(InputStream, current_string);


    if (*KernelVersion==-1)
    {
        *KernelVersion= getIntValue(current_string);
    }
    std::getline(InputStream, current_string);
    if (*KernelVersionSTD==-1)
    {
        *KernelVersionSTD= getIntValue(current_string);
    }
    std::getline(InputStream, current_string);
    *DiscBottom= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *ReadInBottom= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *PartialDryTreatment= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *FrictionTerms= getIntValue(current_string);
    std::getline(InputStream, current_string);
    *ConvertToKM= getIntValue(current_string);
}


void ReadCartesianData(const int fixedDomain,const int fixedDisc, dfloat *xL,dfloat *xR,dfloat *yL,dfloat *yR, int *NelemX, int *NelemY,bool *PeriodicBD_X,bool *PeriodicBD_Y)
{

//,dfloat T, dfloat g_const


    std::ifstream InputStream;
    string filename="Cartesian.dat";
    InputStream.open(filename.c_str());

    if (!InputStream)
    {
        std::string error_message("ERROR: Input file not found: ");
        error_message += filename;
        throw std::invalid_argument(error_message);
    }


    std::string current_string;

    if (!fixedDomain)
    {
        std::getline(InputStream, current_string);
        *xL= getDfloatValue(current_string);
        std::getline(InputStream, current_string);
        *xR= getDfloatValue(current_string);
        std::getline(InputStream, current_string);
        *yL= getDfloatValue(current_string);
        std::getline(InputStream, current_string);
        *yR= getDfloatValue(current_string);
    }
    else
    {
        std::getline(InputStream, current_string);
        std::getline(InputStream, current_string);
        std::getline(InputStream, current_string);
        std::getline(InputStream, current_string);
    }
    if (!fixedDisc)
    {
        std::getline(InputStream, current_string);
        if (*NelemX==0)
        {
            *NelemX= getIntValue(current_string);
        }

        std::getline(InputStream, current_string);
        if (*NelemY==0)
        {
            *NelemY= getIntValue(current_string);
        }

        std::getline(InputStream, current_string);
        *PeriodicBD_X= getBoolValue(current_string);
        std::getline(InputStream, current_string);
        *PeriodicBD_Y= getBoolValue(current_string);
    }
}




void WriteFullMesh(const int NumNodes, const dfloat *x,const dfloat *y)
{

//,dfloat T, dfloat g_const

    const dfloat earth_radius= 6.378e6;
    std::ofstream OutputStream;
    string filename="FullMesh.txt";
    OutputStream.open(filename.c_str());
    for(int inode=0; inode<NumNodes; ++inode)
    {


        dfloat x_degree= (x[inode]/earth_radius)*(180./M_PI);
        dfloat y_degree= (atan(sinh(y[inode]/earth_radius))) * (180./M_PI);
        OutputStream << x_degree << " " << y_degree << "\n";

    }
    OutputStream.close();
}


void ReadFullMesh(const int NumNodes, dfloat *b, dfloat * h_0)
{

//,dfloat T, dfloat g_const
    std::ifstream InputStream;
    string filename="bottomfiles/FullBottomCartesian.txt";
    InputStream.open(filename.c_str());

    if (!InputStream)
    {
        std::string error_message("ERROR: Bottom Topography file not found: ");
        error_message += filename;
        throw std::invalid_argument(error_message);
    }

    std::string current_string;
    std::stringstream current_line(current_string);

    dfloat b_min=0.0;



    for (unsigned i = 0; i < NumNodes; ++i)
    {
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        dfloat dummyA;
        dfloat dummyB;
        dfloat b_input;
        if (!(current_line >> dummyA >> dummyB >> b_input))
        {
            std::string error_message("ERROR: Cant read in Nodes! ");
            error_message += filename;
            throw std::invalid_argument(error_message);
        }
        b[i] = b_input;
        b_min = min(b_min, b_input);
    }
    InputStream.close();
    *h_0 = - floor(b_min);
    for (unsigned i = 0; i < NumNodes; ++i)
    {
        b[i] = *h_0 + b[i];
    }
}
