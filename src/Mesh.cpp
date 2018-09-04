#include "Mesh.h"


Mesh::Mesh(const dfloat * fm_x_GL,const int int_ngl, const int intNele, const int ReadInBottom)
{
//    x_GL=fm_x_GL;
    ngl=int_ngl;
    ngl2 = ngl*ngl;
//    x_GL.resize(ngl);
    x_GL= (dfloat*) calloc(ngl,sizeof(dfloat));
    NelemX=intNele;
    NelemY=intNele;
    ReadBottom = ReadInBottom;

    for (int i=0; i<ngl; i++)
    {

        x_GL[i] = fm_x_GL[i];

    }
//cout << "LGL Nodes: ";
//for(int i = 0; i < ngl; ++i){
//cout << " " << x_GL(i+1) << " " ;
//}
//cout << " \n";
}
Mesh::~Mesh()
{
    //dtor
}


void Mesh::InitMesh(const string meshFile, const bool Cartesian, const int Testcase)
{
    // NelemX=0;
    // NelemY=0;
    if (Cartesian)
    {
        dfloat xR,xL,yR,yL;
        bool PeriodicBD_X,PeriodicBD_Y;
        int fixedDomain,fixedDisc;


        InitDomain(Testcase, &fixedDomain,&fixedDisc,&NelemX, &NelemY,&PeriodicBD_X, &PeriodicBD_Y, &xL, &xR, &yL, &yR);
        ReadCartesianData(fixedDomain,fixedDisc,&xL,&xR,&yL,&yR,&NelemX,&NelemY,&PeriodicBD_X,&PeriodicBD_Y);
        GenerateMesh(xL,xR,yL,yR,PeriodicBD_X,PeriodicBD_Y);
    }
    else
    {
        ReadMesh(meshFile);
        cout << "Mesh read in.\n";
    }




//    ElemEdgeMasterSlave.resize(4,m_num_elements);
//    ElemEdgeOrientation.resize(4,m_num_elements);

    ElemEdgeMasterSlave = (int*) calloc(m_num_elements*4,sizeof(int));
    ElemEdgeOrientation = (int*) calloc(m_num_elements*4,sizeof(int));

    for(int ie=0; ie<m_num_elements; ++ie)
    {
        for (int is=0; is<4; is++)
        {
            int id = ie*4+is;
            int ifa  = ElementToEdge[id];

            if (EdgeInfo[(ifa-1)*7+2]==ie)
            {
                //this is the left element to this edge!
                ElemEdgeMasterSlave[id]=+1;
                ElemEdgeOrientation[id]=1;   //order is never reversed for left element! as it is master

            }
            else
            {
                ElemEdgeMasterSlave[id]=-1;
                ElemEdgeOrientation[id]=EdgeInfo[(ifa-1)*7+6];
            }

        }
    }





    // SPECIAL INTERIOR BOUNDARIES FOR EXAMPLE CURVED DAM BREAK




    if(Testcase ==43 || Testcase == 46 || Testcase == 47)    //PARTIAL CURVED DAM BREAK
    {
        cout << "Doing extra boundaries for Testcase 43, partial curved dam...\n";
        int ExtraEdges=36;
        int ghostCounter = 0;
//        int EdgeInfo_TMP(7,m_num_edges);

        int * EdgeInfo_TMP  = (int*) calloc(m_num_edges*7,sizeof(int));

        for(int is=0; is<m_num_edges; ++is)
        {
            for (int i=0; i<7; i++)
            {
                int id = is*7+i;
                EdgeInfo_TMP[id] = EdgeInfo[id];
            }
        }

        free(EdgeInfo);
        EdgeInfo = (int*) calloc((m_num_edges+ExtraEdges)*7,sizeof(int));
//        EdgeInfo.resize(7,m_num_edges+ExtraEdges);

        for(int is=0; is<m_num_edges; ++is)
        {
            for (int i=0; i<7; i++)
            {
                int id = is*7+i;
                EdgeInfo[id] = EdgeInfo_TMP[id];
            }
        }
//
//
        for(int is=0; is<m_num_edges; ++is)
        {
            int id = is*7;

            if( (EdgeInfo[id+2] % 20 == 2) && (EdgeInfo[id+2] % 40 != 2) && (EdgeInfo[id+3] % 20 == 3)&& (EdgeInfo[id+3] % 40 != 3))
            {
                if ( (is!=1543) && (is!=1624) && (is != 1705) && (is!=1786) )   //&& (is !=1868)
                {


                    int ghostID = (m_num_edges+ghostCounter)*7;
                    ghostCounter = ghostCounter+1;
                    EdgeInfo[ghostID]= EdgeInfo[id];
                    EdgeInfo[ghostID+1]= EdgeInfo[id+1];
                    EdgeInfo[ghostID+2]= EdgeInfo[id+3];
                    EdgeInfo[ghostID+3]= -1;
                    EdgeInfo[ghostID+4]= EdgeInfo[id+5];
                    EdgeInfo[ghostID+5]= -1;
                    EdgeInfo[ghostID+6]= EdgeInfo[id+6];

                    EdgeInfo[id+3] = -1 ;
                    EdgeInfo[id+5]=-1;
                }
            }
        }
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        ghostCounter = 1;
        for(int ie=24; ie<=1584; ie=ie+40)
        {
            if ( (ie!=744) && (ie!=784) && (ie != 824) && (ie!=864) )   //&& (ie !=903)
            {
                int id = (ie-1)*4 + 3;
                ElementToEdge[id]   = m_num_edges + ghostCounter;
                ElemEdgeMasterSlave[id]=+1;
                ElemEdgeOrientation[id]=1;   //order is never reversed for left element! as it is master

                ghostCounter = ghostCounter+1;
            }
        }
//
//
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        m_num_edges=m_num_edges+ExtraEdges;
        cout << "... finished\n";
    }

    if(Testcase ==45)    //PARTIAL CURVED DAM BREAK
    {
        int ExtraEdges = 40;
        int ghostCounter = 0;
        int * EdgeInfo_TMP  = (int*) calloc(m_num_edges*7,sizeof(int));

        for(int is=0; is<m_num_edges; ++is)
        {
            for (int i=0; i<7; i++)
            {
                int id = is*7+i;
                EdgeInfo_TMP[id] = EdgeInfo[id];
            }
        }


        free(EdgeInfo);
        EdgeInfo = (int*) calloc((m_num_edges+ExtraEdges)*7,sizeof(int));
//        EdgeInfo.resize(7,m_num_edges+ExtraEdges);

        for(int is=0; is<m_num_edges; ++is)
        {
            for (int i=0; i<7; i++)
            {
                int id = is*7+i;
                EdgeInfo[id] = EdgeInfo_TMP[id];
            }
        }
//
//
//
        for(int is=0; is<m_num_edges; ++is)
        {
            int id = is*7;
            if( (EdgeInfo[id+2] % 20 == 2) && (EdgeInfo[id+2] % 40 != 2) && (EdgeInfo[id+3] % 20 == 3)&& (EdgeInfo[id+3] % 40 != 3))
            {


                int ghostID = (m_num_edges+ghostCounter)*7;
                ghostCounter = ghostCounter+1;
                EdgeInfo[ghostID]= EdgeInfo[id];
                EdgeInfo[ghostID+1]= EdgeInfo[id+1];
                EdgeInfo[ghostID+2]= EdgeInfo[id+3];
                EdgeInfo[ghostID+3]= -1;
                EdgeInfo[ghostID+4]= EdgeInfo[id+5];
                EdgeInfo[ghostID+5]= -1;
                EdgeInfo[ghostID+6]= EdgeInfo[id+6];

                EdgeInfo[id+3] = -1 ;
                EdgeInfo[id+5]=-1;

            }
        }
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        ghostCounter = 1;
        for(int ie=24; ie<=1584; ie=ie+40)
        {
            int id = (ie-1)*4 + 3;
            ElementToEdge[id]   = m_num_edges + ghostCounter;
            ElemEdgeMasterSlave[id]=+1;
            ElemEdgeOrientation[id]=1;   //order is never reversed for left element! as it is master

            ghostCounter = ghostCounter+1;

        }
//
//
//        cout << "ghostcounter is " << ghostCounter <<"\n";
        m_num_edges=m_num_edges+ExtraEdges;
    }









//    if(Testcase ==20)    //PARTIAL DAM BREAK, VARIABLE NelemX,NelemY
//    {
//        int ExtraEdges= int(0.9*NelemX); // the dam is supposed to have length of 1, so 1/10 of the elements
//        int CutEdges = int(0.1*NelemX);
//        int ghostCounter = 0;
//
//
//        int * EdgeInfo_TMP  = (int*) calloc(m_num_edges*7,sizeof(int));
//
//        for(int is=0; is<m_num_edges; ++is)
//        {
//            for (int i=0; i<7; i++)
//            {
//                int id = is*7+i;
//                EdgeInfo_TMP[id] = EdgeInfo[id];
//            }
//        }
//
//
//        free(EdgeInfo);
//        EdgeInfo = (int*) calloc((m_num_edges+ExtraEdges)*7,sizeof(int));
//
//        for(int is=0; is<m_num_edges; ++is)
//        {
//            for (int i=0; i<7; i++)
//            {
//                int id = is*7+i;
//                EdgeInfo[id] = EdgeInfo_TMP[id];
//            }
//        }

//
//        int startIndex = (NelemX+1)*NelemY + NelemY/2;
//        int endIndex = m_num_edges - NelemY/2 ;
//        int increment =NelemX+1;
//
//        int MidIndex = startIndex + increment*(ExtraEdges/2-1);
//        int startIndex2 = MidIndex + (CutEdges+1) * increment-1;
//

//        for(int is=startIndex; is<MidIndex; is+=increment)
//        {
//            int id = is*7;
//            int ghostID = (m_num_edges+ghostCounter)*7;
//            ghostCounter = ghostCounter+1;
//
//            EdgeInfo[ghostID]= EdgeInfo[id];
//            EdgeInfo[ghostID+1]= EdgeInfo[id+1];
//            EdgeInfo[ghostID+2]= EdgeInfo[id+3];
//            EdgeInfo[ghostID+3]= -1;
//            EdgeInfo[ghostID+4]= EdgeInfo[id+5];
//            EdgeInfo[ghostID+5]= -1;
//            EdgeInfo[ghostID+6]= EdgeInfo[id+6];
//
//            EdgeInfo[id+3] = -1 ;
//            EdgeInfo[id+5]=-1;
//
//
//
//        }
//
//        for(int is=startIndex2; is<endIndex; is+=increment)
//        {
//            int id = is*7;
//            int ghostID = (m_num_edges+ghostCounter)*7;
//            ghostCounter = ghostCounter+1;
//
//            EdgeInfo[ghostID]= EdgeInfo[id];
//            EdgeInfo[ghostID+1]= EdgeInfo[id+1];
//            EdgeInfo[ghostID+2]= EdgeInfo[id+3];
//            EdgeInfo[ghostID+3]= -1;
//            EdgeInfo[ghostID+4]= EdgeInfo[id+5];
//            EdgeInfo[ghostID+5]= -1;
//            EdgeInfo[ghostID+6]= EdgeInfo[id+6];
//
//            EdgeInfo[id+3] = -1 ;
//            EdgeInfo[id+5]=-1;
//
//        }
//        ghostCounter = 1;
//
//        int eleStartIndex = NelemX/2+1;
//        int eleEndIndex = m_num_elements - NelemY/2 +1;
//        int eleIncrement = NelemX;
//        int eleMidIndex = eleStartIndex + eleIncrement*(ExtraEdges/2-1);
//        int eleStartIndex2 = eleMidIndex + (CutEdges+1) * eleIncrement;
//
//
//
//
//
//        for(int ie=eleStartIndex; ie<=eleMidIndex; ie=ie+eleIncrement)
//        {
//            int id = (ie-1)*4 + 3;
//            ElementToEdge[id]   = m_num_edges + ghostCounter;
//            ElemEdgeMasterSlave[id]=+1;
//            ElemEdgeOrientation[id]=1;   //order is never reversed for left element! as it is master
//            ghostCounter = ghostCounter+1;
//        }
//        for(int ie=eleStartIndex2; ie<=eleEndIndex; ie=ie+eleIncrement)
//        {
//            int id = (ie-1)*4 + 3;
//            ElementToEdge[id]   = m_num_edges + ghostCounter;
//            ElemEdgeMasterSlave[id]=+1;
//            ElemEdgeOrientation[id]=1;   //order is never reversed for left element! as it is master
//            ghostCounter = ghostCounter+1;
//        }
//
//
//
//
//        m_num_edges=m_num_edges+ExtraEdges;
//
//    }




    NormalsX = (dfloat*) calloc(m_num_edges*ngl,sizeof(dfloat));
    NormalsY = (dfloat*) calloc(m_num_edges*ngl,sizeof(dfloat));
    Scal = (dfloat*) calloc(m_num_edges*ngl,sizeof(dfloat));
    for (int is = 0; is<m_num_edges; is++)
    {
        for (int i =0; i<ngl; i++)
        {
            int id = is*ngl + i;
            int idLoc = EdgeInfo[is*7+2]*4*ngl + EdgeInfo[is*7+4]*ngl  +i;
            NormalsX[id]	=   nx_global[idLoc];
            NormalsY[id]	=   ny_global[idLoc];
            Scal[id] 	    =   scal_global[idLoc];
            if (EdgeInfo[is*7+3] != -1)
            {
                if (EdgeInfo[is*7+6] != 1)
                {
                    int idright = EdgeInfo[is*7+3]*4*ngl + EdgeInfo[is*7+5]*ngl  +ngl-1 -i;
//		cout << " ID " << is  << " NODE : " << i << " Normals left: ( " << nx_global[idLoc] << ", " << ny_global[idLoc] << ", " << scal_global[idLoc]  <<" )       " ;
//		cout << "Normals right: ( " << nx_global[idright] << ", " << ny_global[idright] << ", " << scal_global[idright]  <<"   " ;
//		cout << "Difference : ( " << nx_global[idright]+nx_global[idLoc] << ", " << ny_global[idLoc]+ny_global[idright] << ", " << scal_global[idLoc]-scal_global[idright]   <<" )  \n " ;

                }
            }
        }

    }
    cout << "Mesh Reading/Generating  completed.\n";
}






void Mesh::GenerateMesh(const dfloat xL,const dfloat xR,const dfloat yL,const dfloat yR,const bool PeriodicBD_X,const bool PeriodicBD_Y)
{

    int ngl2=ngl*ngl;

    cout <<"NelemX: "<<NelemX<< "  NelemY: "<<NelemY <<".\n";
    m_num_elements=NelemX*NelemY;



    m_num_edges=(NelemX+1)*NelemY + (NelemY+1)*NelemX;

    if (PeriodicBD_X)
    {
        m_num_edges = m_num_edges -NelemX;
    }
    if (PeriodicBD_Y)
    {
        m_num_edges = m_num_edges -NelemY;

    }
    dfloat deltaX=(xR-xL)/NelemX;
    dfloat deltaY=(yR-yL)/NelemY;


    dfloat nx,ny;

    x_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    y_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    xXi_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    xEta_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    yXi_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    yEta_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    J_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    nx_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));
    ny_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));
    scal_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));

    EdgeInfo = (int*) calloc(m_num_edges*7,sizeof(int));

    //store global edge number for element local sides 1-4
    ElementToEdge = (int*) calloc(m_num_elements*4,sizeof(int));
//    ElementToEdge.resize(4,m_num_elements);


    dfloat x_xi=0.5*deltaX;
    dfloat y_eta=0.5*deltaY;
    dfloat Jac=x_xi*y_eta;

    //generate mesh
    for(int ieY=0; ieY<NelemY; ++ieY)
    {
        for(int ieX=0; ieX<NelemX; ++ieX)
        {
            for(int j=0; j<ngl; ++j)
            {
                for(int i=0; i<ngl; ++i)
                {
                    int id = (ieY*NelemX+ieX)*ngl2   +j*ngl+i;
                    x_global[id] = xL + (x_GL[i]+1.0)/2.0 * deltaX + ieX* deltaX;
                    y_global[id] = yL + (x_GL[j]+1.0)/2.0 * deltaY + ieY* deltaY;
//                cout <<"ID: " <<id <<" (x,y) = ("<<x_global(id)<<" , "<<y_global(id)<<")\n";
                }
            }
        }
    }



    for (unsigned ie = 0; ie < m_num_elements; ++ie)
    {
        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {

                int id = ie*ngl2   +j*ngl+i;

                xXi_global[id] = x_xi;
                xEta_global[id] = 0.0;
                yXi_global[id] = 0.0;
                yEta_global[id] = y_eta;
                J_global[id] = Jac;
            }
        }


    }

    int jLoopBound = NelemY;

    if (PeriodicBD_X)
    {
        jLoopBound = jLoopBound - 1;
    }
    //eta==const edges first
    for (int j=0; j<=jLoopBound; j++)
    {
        for (int i=0; i<NelemX; i++)
        {
            int is = (j*NelemX+i)+1;

//            cout << "EDGE: " << is <<"\n";
            if (j==0)
            {
                // WE SWAP THE USUAL CONVENTION OF NUMBERING HERE, SUCH THAT THE LEFT ELEMENT ON EACH EDGE IS NOT -1

                EdgeInfo[(is-1)*7+2] = i; //left element
                EdgeInfo[(is-1)*7+4] =0;  // side left element
                EdgeInfo[(is-1)*7+3] =-1; // right element . HERE: outer boundary
                EdgeInfo[(is-1)*7+5] =-1; // side right element . HERE: outer boundary


//                for(int iloc=1;iloc<=ngl;++iloc){
//                    int id = EdgeInfo(5,is)*ngl  +iloc;
//                    ny_global(id,EdgeInfo(3,is)+1) = -1.0;
//                    nx_global(id,EdgeInfo(3,is)+1) = 0.0;
//                    scal_global(id,EdgeInfo(3,is)+1) = x_xi;
//                }
                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    // determine index by EleNumber*StorageAmountPerElement + SideinElement*StoragePerSide + IndexOnSide
                    int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                    ny_global[id] = -1.0;
                    nx_global[id] = 0.0;
                    scal_global[id] = x_xi;
                }





                if (PeriodicBD_X)
                {
                    EdgeInfo[(is-1)*7+3] = (NelemY-1)*NelemX+i;
                    EdgeInfo[(is-1)*7+5] =2;

                    for(int iloc=0; iloc<ngl; ++iloc)
                    {
                        int id = EdgeInfo[(is-1)*7+3]*4*ngl + EdgeInfo[(is-1)*7+5]*ngl  +iloc;
//                        int id = EdgeInfo(6,is)*ngl  +iloc;
                        ny_global[id] = 1.0;
                        nx_global[id] = 0.0;
                        scal_global[id] = x_xi;
                    }
                }




            }
            else if(j==NelemY)
            {


                if (!PeriodicBD_X)
                {

                    EdgeInfo[(is-1)*7+2] =(NelemY-1)*NelemX+i;
                    EdgeInfo[(is-1)*7+4] =2;
                    EdgeInfo[(is-1)*7+3] =-1;
                    EdgeInfo[(is-1)*7+5] =-1;
                    for(int iloc=0; iloc<ngl; ++iloc)
                    {
                        int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                        ny_global[id] = 1.0;
                        nx_global[id] = 0.0;
                        scal_global[id] = x_xi;
//                        cout << "nx_global (" << id << " , " <<EdgeInfo(3,is)+1 <<" ) : " <<nx_global(id,EdgeInfo(3,is)+1)  <<"\n";
                    }

                }




            }
            else
            {
                EdgeInfo[(is-1)*7+2] =(j-1)*NelemX+i;
                EdgeInfo[(is-1)*7+3]= j   *NelemX+i;
                EdgeInfo[(is-1)*7+4] =2;
                EdgeInfo[(is-1)*7+5] =0;
                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                    ny_global[id] = 1.0;
                    nx_global[id] = 0.0;
                    scal_global[id] = x_xi;
                }
                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    int id = EdgeInfo[(is-1)*7+3]*4*ngl + EdgeInfo[(is-1)*7+5]*ngl  +iloc;
                    ny_global[id] = -1.0;
                    nx_global[id] = 0.0;
                    scal_global[id]= x_xi;
                }

            };
            EdgeInfo[(is-1)*7]=0;
            EdgeInfo[(is-1)*7+1]=0;

            EdgeInfo[(is-1)*7+6] =1;

            int id = EdgeInfo[(is-1)*7+2]*4 + EdgeInfo[(is-1)*7+4];
            ElementToEdge[id]=is;
//            ElementToEdge(EdgeInfo(5,is)+1,EdgeInfo(3,is)+1)=is;


            if ( EdgeInfo[(is-1)*7+3] != -1)
            {
//                cout <<"Element: " <<EdgeInfo(4,is)+1 <<"   side 1 is edge  "<<is <<"\n" ;
                int id = EdgeInfo[(is-1)*7+3]*4 + EdgeInfo[(is-1)*7+5];
                ElementToEdge[id]=is;

            }

        }
    }

    int iLoopBound = NelemX;

    if (PeriodicBD_Y)
    {
        iLoopBound = iLoopBound - 1;
    }
    //xi==const edges now
    for (int j=0; j<NelemY; j++)
    {
        for (int i=0; i<=iLoopBound; i++)
        {
            int is = ((jLoopBound+1)*NelemX) + (j*(iLoopBound+1)+i)+1;        //12 + ....       2 +2
//            cout << "EDGE xi: " << is <<"\n";

            if (i==0)
            {
                // WE SWAP THE USUAL CONVENTION OF NUMBERING HERE, SUCH THAT THE LEFT ELEMENT ON EACH EDGE IS NOT -1



                EdgeInfo[(is-1)*7+2]  = j*NelemX;
                EdgeInfo[(is-1)*7+4] =3;
                EdgeInfo[(is-1)*7+3]  = -1;
                EdgeInfo[(is-1)*7+5] =-1;

                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                    ny_global[id] = 0.0;
                    nx_global[id] = -1.0;
                    scal_global[id] = y_eta;
                }

                if (PeriodicBD_Y)
                {

                    EdgeInfo[(is-1)*7+3]  =  (j+1)*NelemX-1;
                    EdgeInfo[(is-1)*7+5]=1;
                    for(int iloc=0; iloc<ngl; ++iloc)
                    {
                        int id = EdgeInfo[(is-1)*7+3]*4*ngl + EdgeInfo[(is-1)*7+5]*ngl  +iloc;
                        ny_global[id] = 0.0;
                        nx_global[id] = 1.0;
                        scal_global[id] = y_eta;
                    }
                }




            }
            else if(i==NelemX)
            {


                // these edges dont exist if Periodic Y boundaries are chosen
                if (!PeriodicBD_Y)
                {
                    EdgeInfo[(is-1)*7+2]  =  (j+1)*NelemX-1;
                    EdgeInfo[(is-1)*7+4] =1;
                    EdgeInfo[(is-1)*7+3]  = -1;
                    EdgeInfo[(is-1)*7+5] =-1;
                    for(int iloc=0; iloc<ngl; ++iloc)
                    {
                        int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                        ny_global[id]= 0.0;
                        nx_global[id] = 1.0;
                        scal_global[id] = y_eta;
                    }
                }







            }
            else
            {
                EdgeInfo[(is-1)*7+2]  =  j*NelemX+i-1;
                EdgeInfo[(is-1)*7+3]  =  j*NelemX+i;
                EdgeInfo[(is-1)*7+4] =1;
                EdgeInfo[(is-1)*7+5] =3;
                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    int id = EdgeInfo[(is-1)*7+2]*4*ngl + EdgeInfo[(is-1)*7+4]*ngl  +iloc;
                    ny_global[id] = 0.0;
                    nx_global[id] = 1.0;
                    scal_global[id]= y_eta;
                }
                for(int iloc=0; iloc<ngl; ++iloc)
                {
                    int id = EdgeInfo[(is-1)*7+3]*4*ngl + EdgeInfo[(is-1)*7+5]*ngl  +iloc;
                    ny_global[id] = 0.0;
                    nx_global[id] = -1.0;
                    scal_global[id] = y_eta;
                }


            };
            EdgeInfo[(is-1)*7]=0;
            EdgeInfo[(is-1)*7+1]=0;

            EdgeInfo[(is-1)*7+6] =1;

//            cout <<"Element: " <<EdgeInfo(3,is)+1 <<"   side 2 is edge  "<<is <<"\n" ;
            int id = EdgeInfo[(is-1)*7+2]*4 + EdgeInfo[(is-1)*7+4];
            ElementToEdge[id]=is;


            if ( EdgeInfo[(is-1)*7+3] != -1)
            {
//            cout <<"Element: " <<EdgeInfo(4,is)+1 <<"   side 4 is edge  "<<is <<"\n" ;
                int id = EdgeInfo[(is-1)*7+3]*4 + EdgeInfo[(is-1)*7+5];
                ElementToEdge[id]=is;

            }
        }
    }

//    for (unsigned ie = 1; ie <= m_num_elements; ++ie){
//        cout <<"Element: "<<ie<<"\n";
//        for (unsigned is=1;is<=4;is++){
//            cout <<"side "<<is<<" is global edge: "<<ElementToEdge(is,ie)<<"\n";
//            cout <<"global edge "<< ElementToEdge(is,ie) <<" has left element "<< EdgeInfo(3,ElementToEdge(is,ie)) <<" and right element "<< EdgeInfo(4,ElementToEdge(is,ie)) <<"\n";
//        }
//
//    }

//for (int is=1;is<=m_num_edges;is++){
//
//    cout << "Edge Number: "<<is<<"   ";
//    cout <<"Left Element: " <<EdgeInfo(3,is)+1 <<"  Right Element:  "<<EdgeInfo(4,is)+1 <<"    " ;
//    cout <<"Edge Left: " <<EdgeInfo(5,is) <<"  Edge Right:  "<<EdgeInfo(6,is)<<"\n" ;
//    cout <<"nx: " ;
//    for(int iloc=1;iloc<=ngl;++iloc){
//        int id = (EdgeInfo(5,is))*ngl  +iloc;
//        cout << " "<<nx_global(id,EdgeInfo(3,is)+1);
////        cout << "nx_global (" << id << " , " <<EdgeInfo(3,is)+1 <<" ) : " <<nx_global(id,EdgeInfo(3,is)+1)  <<"\n";
//    }
//    cout <<"      ny: " ;
//    for(int iloc=1;iloc<=ngl;++iloc){
//        int id = (EdgeInfo(5,is))*ngl  +iloc;
//        cout << " "<<ny_global(id,EdgeInfo(3,is)+1);
//    }
//    cout <<"\n";
//}



    // STORING EDGE INFORMATION!!
    //1 = start node id of edge
    //2 = end node id of edge
    //3 = id of left element
    //4 = id of right element
    //5 = local edge left element
    //6 = local edge right element
    //7 = is orientation swapped?

}

void Mesh::ReadMesh(const string meshFile)
{

//,dfloat T, dfloat g_const

//    corners.resize(4,2);
//    x_phy.resize(ngl,ngl);
//    y_phy.resize(ngl,ngl);









    std::ifstream InputStream;
    string filename=meshFile;
    InputStream.open(filename.c_str());

    if (!InputStream)
    {
        std::string error_message("ERROR: Mesh file not found: ");
        error_message += filename;
        throw std::invalid_argument(error_message);
    }

    std::string FileFormat;
    std::string current_string;
    std::getline(InputStream, current_string);
    std::stringstream current_line(current_string);


    current_line>>FileFormat;
    cout << "Fileformat: "<<FileFormat<<"\n";

    std::getline(InputStream, current_string);
    current_line.clear();
    current_line.str(current_string);


    if (!(current_line >> m_num_nodes >> m_num_edges >> m_num_elements >> m_order_of_boundary_edges))
    {
        std::string error_message("ERROR: Header incorrect! ");
        error_message += filename;
        throw std::invalid_argument(error_message);
    }

    dfloat x_nodes[m_num_nodes];
    dfloat y_nodes[m_num_nodes];
    dfloat b_nodes[m_num_nodes];
    //initialise chebyshev nodes immediately
    x_cheby= (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
    w_bary= (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));


    for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
    {
        x_cheby[k]= -cos(k * PI / m_order_of_boundary_edges);   // k = m_order_of_boundary_edges/2 => x_cheby(k+1) = -cos(PI/2) = 0
    }
    //calculate barycentric weights for chebyshev nodes
    BarycentricWeights();
//cout << "m_order_of_bd_edges " <<m_order_of_boundary_edges <<"\n";
//    for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
//    {
//        cout << "Bary Weights ("<<k<<") : "<<w_bary[k] <<"\n";
//
//    }
//    for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
//    {
//        cout << "x_cheby ("<<k<<") : "<<x_cheby[k] <<"\n";
//
//    }
//cout << "Continuing \n";


    for (unsigned i = 0; i < m_num_nodes; ++i)
    {
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        if (!(current_line >> x_nodes[i] >> y_nodes[i] >> b_nodes[i]))
        {
            std::string error_message("ERROR: Cant read in Nodes! ");
            error_message += filename;
            throw std::invalid_argument(error_message);
        }
    }


    // STORING EDGE INFORMATION!!
    //1 = start node id of edge
    //2 = end node id of edge
    //3 = id of left element
    //4 = id of right element
    //5 = local edge left element
    //6 = local edge right element
    //7 = is orientation swapped?
//    int EdgeData[7*m_num_edges];


    EdgeInfo = (int*) calloc(m_num_edges*7,sizeof(int));
    //store global edge number for element local sides 1-4
    ElementToEdge = (int*) calloc(m_num_elements*4,sizeof(int));
//    ElementToEdge.resize(4,m_num_elements);

    for (int is = 0; is < m_num_edges; ++is)
    {
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        int start_node_id, end_node_id, element_id_on_left, element_id_on_right, side_of_left_element;
        int side_of_right_element;
        if (!(current_line >> start_node_id >> end_node_id
                >> element_id_on_left >> element_id_on_right
                >> side_of_left_element >> side_of_right_element))
        {
            std::string error_message("ERROR: Cant read in Edge Stuff! ");
            error_message += filename;
            throw std::invalid_argument(error_message);
        }
        int idEdgeInfo = is*7;
        EdgeInfo[idEdgeInfo]=start_node_id-1;
        EdgeInfo[idEdgeInfo+1]= end_node_id-1;
        EdgeInfo[idEdgeInfo+2]= element_id_on_left-1;
        EdgeInfo[idEdgeInfo+3]= element_id_on_right-1;
        EdgeInfo[idEdgeInfo+4]= side_of_left_element-1;

        int id = (element_id_on_left-1)*4 + side_of_left_element -1;
        ElementToEdge[id]=is+1;



        // sgn(side_of_right_element)

        int index_direction;
        if (side_of_right_element<0)
        {
            //SIDE IS REVERSED
            index_direction=0;
            EdgeInfo[idEdgeInfo+5]= -side_of_right_element-1;
        }
        else
        {
            index_direction=1;
            EdgeInfo[idEdgeInfo+5]= side_of_right_element-1;
        }
        EdgeInfo[idEdgeInfo+6]= index_direction;

        if (EdgeInfo[idEdgeInfo+3]!=-1)
        {
            if (side_of_right_element<0)
            {
                //SIDE IS REVERSED
                int id = (element_id_on_right-1)*4 - side_of_right_element -1;
                ElementToEdge[id]=is+1;
            }
            else
            {
                int id = (element_id_on_right-1)*4 + side_of_right_element -1;
                ElementToEdge[id]=is+1;
            }
        }


//          ElementToEdge(EdgeInfo(6,is)+1,EdgeInfo(4,is)+1)=is;


    }



//    for (unsigned ie = 1; ie <= m_num_elements; ++ie){
//        cout <<"Element: "<<ie<<"\n";
//        for (unsigned is=1;is<=4;is++){
//            cout <<"side "<<is<<" is global edge: "<<ElementToEdge(is,ie)<<"\n";
//            cout <<"global edge "<< ElementToEdge(is,ie) <<" has left element "<< EdgeInfo(3,ElementToEdge(is,ie)) <<" and right element "<< EdgeInfo(4,ElementToEdge(is,ie)) <<"\n";
//        }
//
//    }


    // STORING ELEMENT DATA INFORMATION
    //1 = Corner node 1     ID!!
    //2 = Corner node 2     ID!!
    //3 = Corner node 3     ID!!
    //4 = Corner node 4     ID!!

    // STORING ELEMENT Boundary Name INFORMATION
    //1 = Boundary Name side 1
    //2 = Boundary Name side 2
    //3 = Boundary Name side 3
    //4 = Boundary Name side 4

//    unsigned eleInfoNo=4;
//    int ElementData[eleInfoNo*m_num_elements];

    string ElementBoundaryNames[4*m_num_elements];



    x_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    y_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    xXi_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    xEta_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    yXi_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    yEta_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    J_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    nx_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));
    ny_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));
    scal_global = (dfloat*) calloc(m_num_elements*4*ngl,sizeof(dfloat));


    // new for read in of bottom data
    if(ReadBottom)
    {
        b_global = (dfloat*) calloc(m_num_elements*ngl2,sizeof(dfloat));
    }

    for (unsigned ie = 0; ie < m_num_elements; ++ie)
    {


        //    Gamma4.resize(m_order_of_boundary_edges+1,2);
        dfloat * Gamma1X = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma1Y = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma2X = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma2Y = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma3X = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma3Y = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma4X = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * Gamma4Y = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        dfloat * cornersX = (dfloat*) calloc(4,sizeof(dfloat));
        dfloat * cornersY = (dfloat*) calloc(4,sizeof(dfloat));

        x_phy= (dfloat*) calloc(ngl2,sizeof(dfloat));
        y_phy= (dfloat*) calloc(ngl2,sizeof(dfloat));
        x_xi= (dfloat*) calloc(ngl2,sizeof(dfloat));
        x_eta= (dfloat*) calloc(ngl2,sizeof(dfloat));
        y_xi= (dfloat*) calloc(ngl2,sizeof(dfloat));
        y_eta= (dfloat*) calloc(ngl2,sizeof(dfloat));
        J= (dfloat*) calloc(ngl2,sizeof(dfloat));
        if(ReadBottom)
        {
            b_phy= (dfloat*) calloc(ngl2,sizeof(dfloat));

            cornersb = (dfloat*) calloc(4,sizeof(dfloat));
        }

        Gamma1b = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        Gamma2b = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        Gamma3b = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));
        Gamma4b = (dfloat*) calloc(m_order_of_boundary_edges+1,sizeof(dfloat));

        x_bndy= (dfloat*) calloc(4*ngl,sizeof(dfloat));
        y_bndy= (dfloat*) calloc(4*ngl,sizeof(dfloat));
        scal= (dfloat*) calloc(4*ngl,sizeof(dfloat));
        nx= (dfloat*) calloc(4*ngl,sizeof(dfloat));
        ny= (dfloat*) calloc(4*ngl,sizeof(dfloat));



        bool curved=false;

        int CornerIDs[4];
        // Get the corner nodes.
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        if (!(current_line >> CornerIDs[0 ]
                >> CornerIDs[1 ]
                >> CornerIDs[2 ]
                >> CornerIDs[3 ]))

        {
            std::string error_message("ERROR: Cant read in Element Corners! ");
            error_message += filename;
            throw std::invalid_argument(error_message);
        }


        // ISM2 format starts counting from one, so we substract one since we count from zero.
        for (unsigned j = 0; j < 4; ++j)
        {
            CornerIDs[j ]-=1;
        }

        // Get the information about whether the boundary curves are given.
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        int is_gamma_given[4];
        if (!(current_line >> is_gamma_given[0] >> is_gamma_given[1] >> is_gamma_given[2] >> is_gamma_given[3]))
        {
            std::string error_message("ERROR: Cant read in Gamma information! ");
            error_message += filename;
            throw std::invalid_argument(error_message);
        }

        // If none of the gammas are given, then we have a fully straight sided quadrilateral.
        // We do not need to set the gammas in this case.
        if (!(is_gamma_given[0] || is_gamma_given[1] || is_gamma_given[2] || is_gamma_given[3]))
        {
            curved=false;
        }//ElementData[eleInfoNo*i + 4 ]=1;    // side is straight sided
        else
        {
            curved=true;
        }//ElementData[eleInfoNo*i + 4 ]=0;    // side is curved

        // Determine if the boundary curve is given.
        // If it is not, then we need to build the boundary curve from the two nodes.
        // If it is, simply copy over the boundary curve information to the intermediate storage class.
        for (unsigned j = 0; j < 4; ++j)
        {
            //initialise the corners arrays or set it to the current corner nodes
            cornersX[j] = x_nodes[CornerIDs[j]];
            cornersY[j] = y_nodes[CornerIDs[j]];
            if(ReadBottom)
            {
                cornersb[j] = b_nodes[CornerIDs[j]];
            }
            //if we dont have a given Gamma curve, we get one by interpolation (chebychev)
            if (is_gamma_given[j] == false)
            {

                unsigned corner_initial;
                unsigned corner_final;

                switch(j)
                {
                case 0:
                    corner_initial= CornerIDs[ 0 ];
                    corner_final= CornerIDs[ 1 ];
                    break;
                case 1:
                    corner_initial= CornerIDs[ 1 ];
                    corner_final= CornerIDs[ 2 ];
                    break;
                case 2:
                    corner_initial= CornerIDs[ 3 ];
                    corner_final= CornerIDs[ 2 ];
                    break;
                case 3:
                    corner_initial= CornerIDs[ 0 ];
                    corner_final= CornerIDs[ 3 ];
                    break;

                }


                dfloat x_initial = x_nodes[corner_initial];
                dfloat x_final = x_nodes[corner_final];

                dfloat y_initial = y_nodes[corner_initial];
                dfloat y_final = y_nodes[corner_final];
                if(ReadBottom)
                {
                    b_initial = b_nodes[corner_initial];
                    b_final = b_nodes[corner_final];
                }
                for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
                {


                    dfloat interpolated_x = ((x_final - x_initial) * x_cheby[k]+ x_initial + x_final)/2.0;
                    dfloat interpolated_y = ((y_final - y_initial) * x_cheby[k] + y_initial + y_final)/2.0;
                    if(ReadBottom)
                    {
                        interpolated_b = ((b_final - b_initial) * x_cheby[k] + b_initial + b_final)/2.0;
                    }

                    switch(j)
                    {
                    case 0:
                        Gamma1X[k]= interpolated_x;
                        Gamma1Y[k]= interpolated_y;
                        if(ReadBottom)
                        {
                            Gamma1b[k] = interpolated_b;
                        }
                        break;
                    case 1:
                        Gamma2X[k]= interpolated_x;
                        Gamma2Y[k]= interpolated_y;
                        if(ReadBottom)
                        {
                            Gamma2b[k] = interpolated_b;
                        }
                        break;
                    case 2:
                        Gamma3X[k]= interpolated_x;
                        Gamma3Y[k]= interpolated_y;
                        if(ReadBottom)
                        {
                            Gamma3b[k] = interpolated_b;
                        }
                        break;
                    case 3:
                        Gamma4X[k]= interpolated_x;
                        Gamma4Y[k]= interpolated_y;
                        if(ReadBottom)
                        {
                            Gamma4b[k] = interpolated_b;
                        }
                        break;

                    }
                }


            }
            else // is_gamma_given is true
            {
//		cout << "ACTUALLY READING IN GAMMA CURVE!\n";
                curved =true;
                for (unsigned k = 0; k <= m_order_of_boundary_edges; ++k)
                {
                    std::getline(InputStream, current_string);
                    current_line.clear();
                    current_line.str(current_string);
                    switch(j)
                    {
                    case 0:
                        if (!(current_line >> Gamma1X[k]>> Gamma1Y[k]>> Gamma1b[k]))
                        {
                            std::string error_message("ERROR: Cant read in given Gamma1 Curve! ");
                            error_message += filename;
                            cout << "Error in Line : "<<k <<"\n";
                            throw std::invalid_argument(error_message);
                        }
                        break;
                    case 1:
                        if (!(current_line >> Gamma2X[k]>> Gamma2Y[k]>> Gamma2b[k]))
                        {
                            std::string error_message("ERROR: Cant read in given Gamma2 Curve! ");
                            error_message += filename;
                            cout << "Error in Line : "<<k <<"\n";
                            throw std::invalid_argument(error_message);
                        }
                        break;
                    case 2:
                        if (!(current_line >> Gamma3X[k]>> Gamma3Y[k]>> Gamma3b[k]))
                        {
                            std::string error_message("ERROR: Cant read in given Gamma3 Curve! ");
                            error_message += filename;
                            cout << "Error in Line : "<<k <<"\n";
                            throw std::invalid_argument(error_message);
                        }
                        break;
                    case 3:
                        if (!(current_line >> Gamma4X[k]>> Gamma4Y[k]>> Gamma4b[k]))
                        {
                            std::string error_message("ERROR: Cant read in given Gamma4 Curve! ");
                            error_message += filename;
                            cout << "Error in Line : "<<k <<"\n";
                            throw std::invalid_argument(error_message);
                        }
                        break;

                    }


                }
            }
        }


        // Extract the boundary names
        std::getline(InputStream, current_string);
        current_line.clear();
        current_line.str(current_string);
        for (unsigned j = 0; j < 4; ++j)
        {
            if (!(current_line >> ElementBoundaryNames[4*ie +j]))
            {
                std::string error_message("ERROR: Cant read in given Boundary Names! ");
                error_message += filename;
                throw std::invalid_argument(error_message);
            }
        }
//	curved =true;
        ConstructMappedGeometry(cornersX,cornersY,cornersb,
                                Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,
                                Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,
                                Gamma1b,Gamma2b,Gamma3b,Gamma4b,
                                curved);


        for(int j=0; j<ngl; ++j)
        {
            for(int i=0; i<ngl; ++i)
            {
                int locID = j*ngl+i;
                int id = ie*ngl2   +locID;

                x_global[id] = x_phy[locID];
                y_global[id] = y_phy[locID];
                xXi_global[id] = x_xi[locID];
                xEta_global[id] = x_eta[locID];
                yXi_global[id] = y_xi[locID];
                yEta_global[id] = y_eta[locID];
                J_global[id] = J[locID];
                if(ReadBottom)
                {
                    b_global[id] = b_phy[locID];
                }
            }
        }

        for(int is=0; is<4; ++is)
        {
            for(int i=0; i<ngl; ++i)
            {

                int id = is*ngl  +i;
                int eleLocID = ie * ngl * 4 + id;

                nx_global[eleLocID] = nx[id];
                ny_global[eleLocID] = ny[id];
                scal_global[eleLocID] = scal[id];
            }
        }



        free(Gamma1X);
        free(Gamma1Y);
        free(Gamma2X);
        free(Gamma2Y);
        free(Gamma3X);
        free(Gamma3Y);
        free(Gamma4X);
        free(Gamma4Y);
        free(cornersX);
        free(cornersY);
        free(x_phy);
        free(y_phy);
        free(x_xi);
        free(x_eta);
        free(y_xi);
        free(y_eta);
        free(J);
        free(x_bndy);
        free(y_bndy);
        free(scal);
        free(nx);
        free(ny);

        if (ReadBottom)
        {
            free(b_phy);
        }


    }//ele loop

    free(w_bary);
    free(x_cheby);


}




void Mesh :: TransfiniteQuadMap(const dfloat * Gamma1X,const dfloat * Gamma1Y,const dfloat * Gamma2X,const dfloat * Gamma2Y,
                                const dfloat * Gamma3X,const dfloat * Gamma3Y,const dfloat * Gamma4X,const dfloat * Gamma4Y,
                                const dfloat psi,const dfloat eta,dfloat *x_out,dfloat *y_out)
{
    dfloat x1_ref, y1_ref,x1_cp, y1_cp;
    dfloat x2_ref, y2_ref,x2_cp, y2_cp;
    dfloat x3_ref, y3_ref,x3_cp, y3_cp;
    dfloat x4_ref, y4_ref,x4_cp, y4_cp;

    EvaluateAt(Gamma1X,Gamma1Y,-1.0,&x1_ref, &y1_ref);
    EvaluateAt(Gamma1X,Gamma1Y, 1.0,&x2_ref, &y2_ref);
    EvaluateAt(Gamma3X,Gamma3Y, 1.0,&x3_ref, &y3_ref);
    EvaluateAt(Gamma3X,Gamma3Y,-1.0,&x4_ref, &y4_ref);

    EvaluateAt(Gamma1X,Gamma1Y,psi,&x1_cp, &y1_cp);
    EvaluateAt(Gamma2X,Gamma2Y,eta,&x2_cp, &y2_cp);
    EvaluateAt(Gamma3X,Gamma3Y,psi,&x3_cp, &y3_cp);
    EvaluateAt(Gamma4X,Gamma4Y,eta,&x4_cp, &y4_cp);

    *x_out = 0.5*((1.0-psi)*x4_cp+(1.0+psi)*x2_cp+(1.0-eta)*x1_cp+(1.0+eta)*x3_cp) - 0.25*((1.0-psi)*((1.0-eta)*x1_ref+(1.0+eta)*x4_ref)+(1.0+psi)*((1.0-eta)*x2_ref+(1.0+eta)*x3_ref));

    *y_out = 0.5*((1.0-psi)*y4_cp+(1.0+psi)*y2_cp+(1.0-eta)*y1_cp+(1.0+eta)*y3_cp) - 0.25*((1.0-psi)*((1.0-eta)*y1_ref+(1.0+eta)*y4_ref)  +(1.0+psi)*((1.0-eta)*y2_ref+(1.0+eta)*y3_ref));


}





void Mesh :: QuadMap(const dfloat * cornersX,const dfloat * cornersY,const dfloat psi,const dfloat eta,dfloat *x_out,dfloat *y_out)
{
// Mapping of the reference square to a straight sided quadrilateral
    *x_out = 0.25*(cornersX[0]*(1.0-psi)*(1.0-eta)+cornersX[1]*(1.0+psi)*(1.0-eta)+ cornersX[2]*(1.0+psi)*(1.0+eta)+cornersX[3]*(1.0-psi)*(1.0+eta));
    *y_out = 0.25*(cornersY[0]*(1.0-psi)*(1.0-eta)+cornersY[1]*(1.0+psi)*(1.0-eta)+ cornersY[2]*(1.0+psi)*(1.0+eta)+cornersY[3]*(1.0-psi)*(1.0+eta));



}

void Mesh :: TransfiniteQuadMetrics(const dfloat * Gamma1X,const dfloat * Gamma1Y,const dfloat * Gamma2X,const dfloat * Gamma2Y,
                                    const dfloat * Gamma3X,const dfloat * Gamma3Y,const dfloat * Gamma4X,const dfloat * Gamma4Y,
                                    const dfloat psi,const dfloat eta,dfloat *X_psi,dfloat *X_eta,dfloat *Y_psi,dfloat *Y_eta)
{
// Computation of the metric terms on a curve-bounded quadrilateral

    dfloat * Xref1 = (dfloat*) calloc(4,sizeof(dfloat));
    dfloat * Xcomp1 = (dfloat*) calloc(4,sizeof(dfloat));
    dfloat * Xpcomp1 = (dfloat*) calloc(4,sizeof(dfloat));
    dfloat * Xref2 = (dfloat*) calloc(4,sizeof(dfloat));
    dfloat * Xcomp2 = (dfloat*) calloc(4,sizeof(dfloat));
    dfloat * Xpcomp2 = (dfloat*) calloc(4,sizeof(dfloat));




    EvaluateAt(Gamma1X,Gamma1Y,-1.0,&Xref1[0],&Xref2[0]);
    EvaluateAt(Gamma1X,Gamma1Y, 1.0,&Xref1[1],&Xref2[1]);
    EvaluateAt(Gamma3X,Gamma3Y, 1.0,&Xref1[2],&Xref2[2]);
    EvaluateAt(Gamma3X,Gamma3Y,-1.0,&Xref1[3],&Xref2[3]);



    EvaluateAt(Gamma1X,Gamma1Y,psi,&Xcomp1[0],&Xcomp2[0]);
    EvaluateAt(Gamma2X,Gamma2Y,eta,&Xcomp1[1],&Xcomp2[1]);
    EvaluateAt(Gamma3X,Gamma3Y,psi,&Xcomp1[2],&Xcomp2[2]);
    EvaluateAt(Gamma4X,Gamma4Y,eta,&Xcomp1[3],&Xcomp2[3]);



    DerivativeAt(Gamma1X,Gamma1Y,psi,&Xpcomp1[0],&Xpcomp2[0]);
    DerivativeAt(Gamma2X,Gamma2Y,eta,&Xpcomp1[1],&Xpcomp2[1]);
    DerivativeAt(Gamma3X,Gamma3Y,psi,&Xpcomp1[2],&Xpcomp2[2]);
    DerivativeAt(Gamma4X,Gamma4Y,eta,&Xpcomp1[3],&Xpcomp2[3]);




    *X_psi = 0.5*(Xcomp1[1]-Xcomp1[3]+(1.0-eta)*Xpcomp1[0]+(1.0+eta)*Xpcomp1[2]) -  0.25*((1.0-eta)*(Xref1[1]-Xref1[0])+(1.0+eta)*(Xref1[2]-Xref1[3]));

    *Y_psi = 0.5*(Xcomp2[1]-Xcomp2[3]+(1.0-eta)*Xpcomp2[0]+(1.0+eta)*Xpcomp2[2]) -0.25*((1.0-eta)*(Xref2[1]-Xref2[0])+(1.0+eta)*(Xref2[2]-Xref2[3]));

    *X_eta = 0.5*((1.0-psi)*Xpcomp1[3]+(1.0+psi)*Xpcomp1[1]+Xcomp1[2]-Xcomp1[0]) -  0.25*((1.0-psi)*(Xref1[3]-Xref1[0])+(1.0+psi)*(Xref1[2]-Xref1[1]));

    *Y_eta = 0.5*((1.0-psi)*Xpcomp2[3]+(1.0+psi)*Xpcomp2[1]+Xcomp2[2]-Xcomp2[0]) - 0.25*((1.0-psi)*(Xref2[3]-Xref2[0])+(1.0+psi)*(Xref2[2]-Xref2[1]));

//dfloat DIFFERENCE1 = fabs(*X_psi - 0.25*((1.0-eta)*(Xref1[1]-Xref1[0])+(1.0+eta)*(Xref1[2]-Xref1[3])) );
//dfloat DIFFERENCE2 = fabs(*Y_psi - 0.25*((1.0-eta)*(Xref2[1]-Xref2[0])+(1.0+eta)*(Xref2[2]-Xref2[3])) );
//dfloat DIFFERENCE3 =  fabs(*X_eta -  0.25*((1.0-psi)*(Xref1[3]-Xref1[0])+(1.0+psi)*(Xref1[2]-Xref1[1])) );
//dfloat DIFFERENCE4 = fabs(*Y_eta -  0.25*((1.0-psi)*(Xref2[3]-Xref2[0])+(1.0+psi)*(Xref2[2]-Xref2[1])) );

//if (max(max(DIFFERENCE1, DIFFERENCE2), max(DIFFERENCE3, DIFFERENCE4)) > pow(10.0,-12)){
//	cout << "DIFF1: " << DIFFERENCE1  << " DIFF2: " << DIFFERENCE2  << " DIFF3: " << DIFFERENCE3  << " DIFF4: " << DIFFERENCE4  << "\n";
//}
//    *X_psi =   0.25*((1.0-eta)*(Xref1[1]-Xref1[0])+(1.0+eta)*(Xref1[2]-Xref1[3]));
//    *Y_psi = 0.25*((1.0-eta)*(Xref2[1]-Xref2[0])+(1.0+eta)*(Xref2[2]-Xref2[3]));
//    *X_eta =  0.25*((1.0-psi)*(Xref1[3]-Xref1[0])+(1.0+psi)*(Xref1[2]-Xref1[1]));
//    *Y_eta =  0.25*((1.0-psi)*(Xref2[3]-Xref2[0])+(1.0+psi)*(Xref2[2]-Xref2[1]));


    free(Xref1);
    free(Xcomp1);
    free(Xpcomp1);
    free(Xref2);
    free(Xcomp2);
    free(Xpcomp2);
}


void Mesh :: TransfiniteQuadMapSingle(const dfloat * Gamma1X,const dfloat * Gamma2X,
                                      const dfloat * Gamma3X,const dfloat * Gamma4X,
                                      const dfloat psi,const dfloat eta,dfloat *x_out)
{
    dfloat x1_ref, y1_ref,x1_cp, y1_cp;
    dfloat x2_ref, y2_ref,x2_cp, y2_cp;
    dfloat x3_ref, y3_ref,x3_cp, y3_cp;
    dfloat x4_ref, y4_ref,x4_cp, y4_cp;

    EvaluateAtSingle(Gamma1X,-1.0,&x1_ref);
    EvaluateAtSingle(Gamma1X, 1.0,&x2_ref);
    EvaluateAtSingle(Gamma3X, 1.0,&x3_ref);
    EvaluateAtSingle(Gamma3X,-1.0,&x4_ref);

    EvaluateAtSingle(Gamma1X,psi,&x1_cp);
    EvaluateAtSingle(Gamma2X,eta,&x2_cp);
    EvaluateAtSingle(Gamma3X,psi,&x3_cp);
    EvaluateAtSingle(Gamma4X,eta,&x4_cp);

    *x_out = 0.5*((1.0-psi)*x4_cp+(1.0+psi)*x2_cp+(1.0-eta)*x1_cp+(1.0+eta)*x3_cp) - 0.25*((1.0-psi)*((1.0-eta)*x1_ref+(1.0+eta)*x4_ref)+(1.0+psi)*((1.0-eta)*x2_ref+(1.0+eta)*x3_ref));


}
void Mesh :: QuadMapSingle(const dfloat * cornersX,const dfloat psi,const dfloat eta,dfloat *x_out)
{
// Mapping of the reference square to a straight sided quadrilateral
    *x_out = 0.25*(cornersX[0]*(1.0-psi)*(1.0-eta)+cornersX[1]*(1.0+psi)*(1.0-eta)+ cornersX[2]*(1.0+psi)*(1.0+eta)+cornersX[3]*(1.0-psi)*(1.0+eta));


}




void Mesh :: QuadMapMetrics(const dfloat * cornersX,const dfloat * cornersY,const dfloat xi,const dfloat eta,dfloat *X_psi,dfloat *X_eta,dfloat *Y_psi,dfloat *Y_eta)
{
// Metric terms on a straight sided quadrilateral

//we assume we have the correct corners stored in the corner arrays
    *X_psi = 0.25*((1.0-eta)*(cornersX[1]-cornersX[0])+(1.0+eta)*(cornersX[2]-cornersX[3]));
    *Y_psi = 0.25*((1.0-eta)*(cornersY[1]-cornersY[0])+(1.0+eta)*(cornersY[2]-cornersY[3]));
    *X_eta = 0.25*((1.0-xi)*(cornersX[3]-cornersX[0])+(1.0+xi)*(cornersX[2]-cornersX[1]));
    *Y_eta = 0.25*((1.0-xi)*(cornersY[3]-cornersY[0])+(1.0+xi)*(cornersY[2]-cornersY[1]));

}



void Mesh :: ConstructMappedGeometry(const dfloat * cornersX,const dfloat * cornersY,const dfloat * cornersb,
                                     const dfloat * Gamma1X,const dfloat * Gamma1Y,const dfloat * Gamma2X,const dfloat * Gamma2Y,
                                     const dfloat * Gamma3X,const dfloat * Gamma3Y,const dfloat * Gamma4X,const dfloat * Gamma4Y,
                                     const dfloat * Gamma1b,const dfloat * Gamma2b,const dfloat * Gamma3b,const dfloat * Gamma4b,
                                     const bool Curved)
{
// Constructor of geometry and metric terms for quadrilateral domains




    dfloat xXi,xEta,yXi,yEta,Jtemp;

    dfloat xXi_2[ngl2],xEta_2[ngl2],yXi_2[ngl2],yEta_2[ngl2];
    dfloat x_tmp, y_tmp;

    for (int j = 0; j<ngl; j++)
    {
        for (int i = 0; i<ngl; i++)
        {
            int ij = j*ngl+i;
            if (Curved)
            {
                TransfiniteQuadMap(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],x_GL[j],&x_phy[ij],&y_phy[ij]);

                if(ReadBottom)
                {
                    TransfiniteQuadMapSingle(Gamma1b,Gamma2b,Gamma3b,Gamma4b,x_GL[i],x_GL[j],&b_phy[ij]);
                }
                TransfiniteQuadMetrics(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],x_GL[j],&x_xi[ij],&x_eta[ij],&y_xi[ij],&y_eta[ij]);
                //QuadMapMetrics(cornersX,cornersY,x_GL[i],x_GL[j],&xXi_2[ij],&xEta_2[ij],&yXi_2[ij],&yEta_2[ij]);
                // QuadMapMetrics(cornersX,cornersY,x_GL[i],x_GL[j],&x_xi[ij],&x_eta[ij],&y_xi[ij],&y_eta[ij]);

//dfloat DIFFERENCE1 = fabs(x_xi[ij] -xXi_2[ij]);
//dfloat DIFFERENCE2 = fabs(x_eta[ij] -xEta_2[ij]);
//dfloat DIFFERENCE3 =  fabs(y_xi[ij] -yXi_2[ij]);
//dfloat DIFFERENCE4 = fabs(y_eta[ij] -yEta_2[ij]);

//if (max(max(DIFFERENCE1, DIFFERENCE2), max(DIFFERENCE3, DIFFERENCE4)) > pow(10.0,-12)){
//	cout << "DIFF1: " << DIFFERENCE1  << " DIFF2: " << DIFFERENCE2  << " DIFF3: " << DIFFERENCE3  << " DIFF4: " << DIFFERENCE4  << "\n";
//}

            }
            else
            {
                QuadMap(cornersX,cornersY,x_GL[i],x_GL[j],&x_phy[ij],&y_phy[ij]);
                if(ReadBottom)
                {
                    QuadMapSingle(cornersb,x_GL[i],x_GL[j],&b_phy[ij]);
                }
                QuadMapMetrics(cornersX,cornersY,x_GL[i],x_GL[j],&x_xi[ij],&x_eta[ij],&y_xi[ij],&y_eta[ij]);
            }
            J[ij]= x_xi[ij]*y_eta[ij]-x_eta[ij]*y_xi[ij];
        }
    }





    for (int j = 0; j<ngl; j++)
    {
        int idSide2 = ngl+j;
        int idSide4 = 3*ngl+j;
        if (Curved)
        {
//            TransfiniteQuadMap(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,1.0,x_GL[j],&x_bndy[idSide2],&y_bndy[idSide2]);
//            QuadMapMetrics(cornersX,cornersY,1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
            TransfiniteQuadMetrics(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
        }
        else
        {
//            QuadMap(cornersX,cornersY,1.0,x_GL[j],&x_bndy[idSide2],&y_bndy[idSide2]);
            QuadMapMetrics(cornersX,cornersY,1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
        }
        Jtemp = xXi*yEta-xEta*yXi;
        scal[idSide2]     = sqrt(yEta*yEta + xEta*xEta);
        nx[idSide2] = copysign(1.0,Jtemp)*(yEta/scal[idSide2]);
        ny[idSide2] = copysign(1.0,Jtemp)*(-xEta/scal[idSide2]);
        if (Curved)
        {
//           TransfiniteQuadMap(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,-1.0,x_GL[j],&x_bndy[idSide4],&y_bndy[idSide4]);
//            QuadMapMetrics(cornersX,cornersY,-1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
            TransfiniteQuadMetrics(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,-1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
        }
        else
        {
//            QuadMap(cornersX,cornersY,-1.0,x_GL[j],&x_bndy[idSide4],&y_bndy[idSide4]);
            QuadMapMetrics(cornersX,cornersY,-1.0,x_GL[j],&xXi,&xEta,&yXi,&yEta);
        }
        Jtemp = xXi*yEta-xEta*yXi;
        scal[idSide4]      =  sqrt(yEta*yEta + xEta*xEta);
        nx[idSide4] = -copysign(1.0,Jtemp)*(yEta/scal[idSide4]);
        ny[idSide4] = -copysign(1.0,Jtemp)*(-xEta/scal[idSide4]);
    }
    for (int i = 0; i<ngl; i++)
    {
        int idSide1 = i;
        int idSide3 = 2*ngl+i;
        if (Curved)
        {
//            TransfiniteQuadMap(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],-1.0,&x_bndy[idSide1],&y_bndy[idSide1]);
//            QuadMapMetrics(cornersX,cornersY,x_GL[i],-1.0,&xXi,&xEta,&yXi,&yEta);
            TransfiniteQuadMetrics(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],-1.0,&xXi,&xEta,&yXi,&yEta);
        }
        else
        {
//            QuadMap(cornersX,cornersY,x_GL[i],-1.0,&x_bndy[idSide1],&y_bndy[idSide1]);
            QuadMapMetrics(cornersX,cornersY,x_GL[i],-1.0,&xXi,&xEta,&yXi,&yEta);
        }
        Jtemp = xXi*yEta-xEta*yXi;
        scal[idSide1]      =  sqrt(yXi*yXi + xXi*xXi);
        nx[idSide1]= -copysign(1.0,Jtemp)*(-yXi/scal[idSide1]);
        ny[idSide1] = -copysign(1.0,Jtemp)*(xXi/scal[idSide1]);
        if (Curved)
        {
//            TransfiniteQuadMap(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],1.0,&x_bndy[idSide3],&y_bndy[idSide3]);
//            QuadMapMetrics(cornersX,cornersY,x_GL[i],1.0,&xXi,&xEta,&yXi,&yEta);
            TransfiniteQuadMetrics(Gamma1X,Gamma1Y,Gamma2X,Gamma2Y,Gamma3X,Gamma3Y,Gamma4X,Gamma4Y,x_GL[i],1.0,&xXi,&xEta,&yXi,&yEta);
        }
        else
        {
//            QuadMap(cornersX,cornersY,x_GL[i],1.0,&x_bndy[idSide3],&y_bndy[idSide3]);
            QuadMapMetrics(cornersX,cornersY,x_GL[i],1.0,&xXi,&xEta,&yXi,&yEta);
        }
        Jtemp = xXi*yEta-xEta*yXi;
        scal[idSide3]      = sqrt(yXi*yXi + xXi*xXi);
        nx[idSide3] = copysign(1.0,Jtemp)*(-yXi/scal[idSide3]);
        ny[idSide3] = copysign(1.0,Jtemp)*(xXi/scal[idSide3]);
    }

}























void Mesh :: BarycentricWeights()
{


    for(int i=0; i<=m_order_of_boundary_edges; i++)
    {
        w_bary[i]=1.0;
    };


    for (int j=1; j<=m_order_of_boundary_edges; j++)
    {
        for (int k=0; k<j; k++)
        {
            w_bary[k]=w_bary[k]*(x_cheby[k]-x_cheby[j]);
            w_bary[j]=w_bary[j]*(x_cheby[j]-x_cheby[k]);
        };
    };

    for (int j=0; j<=m_order_of_boundary_edges; j++)
    {
        w_bary[j]=1.0/w_bary[j];
    };
//cout <<"bary weights computed!\n";
};




void Mesh :: EvaluateAt(const dfloat * GammaX,const dfloat * GammaY,const dfloat s,dfloat *x_point,dfloat *y_point)
{
// Evaluates a member of the CurveInterpolant type at a point s

    LagrangeInterpolation(s,GammaX,x_point);
    LagrangeInterpolation(s,GammaY,y_point);
}

void Mesh :: DerivativeAt(const dfloat * GammaX,const dfloat * GammaY,const dfloat s,dfloat *x_point_prime,dfloat *y_point_prime)
{
// Evaluates the derivative of a member of the CurveInterpolant type at a point s


    LagrangeInterpolantDerivative(s,GammaX,x_point_prime);
    LagrangeInterpolantDerivative(s,GammaY,y_point_prime);

}

void Mesh :: EvaluateAtSingle(const dfloat * GammaX,const dfloat s,dfloat *x_point)
{
// Evaluates a member of the CurveInterpolant type at a point s

    LagrangeInterpolation(s,GammaX,x_point);
}

void Mesh :: DerivativeAtSingle(const dfloat * GammaX,const dfloat s,dfloat *x_point_prime)
{
// Evaluates the derivative of a member of the CurveInterpolant type at a point s


    LagrangeInterpolantDerivative(s,GammaX,x_point_prime);

}





void Mesh :: LagrangeInterpolantDerivative(const dfloat xpt,const dfloat* functionvals,dfloat *p_prime)
{

    bool atNode      = false;
    dfloat numerator   = 0.0;
    dfloat denominator=0.0;
    dfloat p,t;
    int k;

    for (int j = 0; j<=m_order_of_boundary_edges; j++)
    {
        if(fabs(xpt-x_cheby[j])<pow(10.0,-14))
        {
            atNode = true;
            p = functionvals[j];
            denominator = -w_bary[j];
            k = j;

        }
    }
    if (atNode)
    {
        for (int j = 0; j<=m_order_of_boundary_edges; j++)
        {
            if(j!=k)
            {
                numerator = numerator + w_bary[j]*(p-functionvals[j])/(xpt-x_cheby[j]);
            }
        }
    }
    else
    {
        denominator = 0.0;
        LagrangeInterpolation(xpt,functionvals,&p);
//	cout << "Interpolated Gamma Curve at " << p << "\n";
        for (int j = 0; j<=m_order_of_boundary_edges; j++)
        {
            t = w_bary[j]/(xpt-x_cheby[j]);
            numerator +=  t*(p-functionvals[j])/(xpt-x_cheby[j]);
            denominator += t;
        }
    }
    *p_prime = numerator/denominator;

}


void Mesh :: LagrangeInterpolation(const dfloat xpt,const dfloat* functionvals,dfloat *output)
{
// Barycentric two formulation of Lagrange interpolant



    dfloat numerator,denominator,t;
    bool var1=false;

    numerator   = 0.0;
    denominator = 0.0;

    for (int j = 0; j<=m_order_of_boundary_edges; j++)
    {
        if(fabs(xpt-x_cheby[j])<pow(10.0,-14))
        {
            var1=true;
            *output = functionvals[j];
        }
    }
    if (!var1)
    {
        for (int j = 0; j<=m_order_of_boundary_edges; j++)
        {
            t = w_bary[j]/(xpt-x_cheby[j]);
            numerator +=  t*functionvals[j];
            denominator +=  t;
        }

        *output = numerator/denominator;
    }

}












void Mesh::InitDomain(const int Testcase, int *fixedDomain,int *fixedDisc, int *NelemX,int *NelemY,bool *PeriodicBD_X, bool*PeriodicBD_Y, dfloat *xL, dfloat *xR, dfloat *yL, dfloat *yR)
{

    switch(Testcase)
    {
    case 1:      // convergence test, periodic
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0 ;
        *xL=0.0;
        *xR=2.0;
        *yL=0.0;
        *yR=2.0;
        break;
    }
    case 7:      // convergence test, periodic	No Botttom
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0 ;
        *xL=0.0;
        *xR=2.0;
        *yL=0.0;
        *yR=2.0;
        break;
    }

    case 2:      // WELL BALANCED (CARTESIAN 20x20)
    {
        *fixedDomain = 0 ;
        *fixedDisc = 1 ;
        *NelemX=20;
        *NelemY=20;
        *PeriodicBD_X=0;
        *PeriodicBD_Y=0;
        break;
    }
    case 3:     // Entropy Glitch
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;
        //*PeriodicBD_X=0;
        //*PeriodicBD_Y=0;
        //*NelemX=40;
        //*NelemY=40;

        *xL=-1.0;
        *xR=1.0;
        *yL=-1.0;
        *yR=1.0;
        break;
    }
    case 9:     // Entropy Glitch	(other way)
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-5.0;
        *xR=5.0;
        *yL=-5.0;
        *yR=5.0;
        break;
    }

    case 4:      // WELL BALANCED (CARTESIAN 20x20)
    {
        *fixedDomain = 0 ;
        *fixedDisc = 1 ;
        *NelemX=4;
        *NelemY=4;
        *PeriodicBD_X=0;
        *PeriodicBD_Y=0;
        break;
    }
    case 5:      // Dry Lake + Lake at Rest
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0 ;
        *xL=0.0;
        *xR=1.0;
        *yL=0.0;
        *yR=1.0;

        break;
    }
    case 20:     // Dam Break
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-5.0;
        *xR=5.0;
        *yL=-5.0;
        *yR=5.0;
        break;
    }
    case 21:     // Dam Break
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-5.0;
        *xR=5.0;
        *yL=-5.0;
        *yR=5.0;
        break;
    }
    case 30:     // Steeper Dam Break To Test Shock Capturing
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-20.0;
        *xR=20.0;
        *yL=-20.0;
        *yR=20.0;
        break;
    }
    case 31:     // Two-dimensional oscillating lake  (Xing_PosPres paper, 6.8)
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-2.0;
        *xR=2.0;
        *yL=-2.0;
        *yR=2.0;


        break;

    }
    case 32:     // Three Mound (4.6)
    {
        *fixedDomain = 1 ;
        *fixedDisc = 1;
        *xL=0.0;
        *xR=75.0;
        *yL=0.0;
        *yR=30.0;
        *NelemX=150;
        *NelemY=100;
        *PeriodicBD_X=0;
        *PeriodicBD_Y=0;


        break;

    }

    case 33:     // Dam Break Three Mound (4.6)
    {

        *fixedDomain = 1 ;
        *fixedDisc = 1;
        *xL=0.0;
        *xR=75.0;
        *yL=0.0;
        *yR=30.0;
        *NelemX=150;
        *NelemY=100;
        *PeriodicBD_X=0;
        *PeriodicBD_Y=0;


        break;

    }


    case 34:      // 2D Solitary Wave Runup and Run-down
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=0.0;
        *xR=25.0;
        *yL=0.0;
        *yR=30.0;


        break;
    }
    case 35:      // 1D Bowl
    {
        *fixedDomain = 1 ;
        *fixedDisc = 0;

        *xL=-3000.0;
        *xR=3000.0;
        *yL=-3000.0;
        *yR=3000.0;


        break;
    }
    case 36:     // Dam Break Three Mound (4.6)
    {

        *fixedDomain = 1 ;
        *fixedDisc = 1;
        *xL=0.0;
        *xR=75.0;
        *yL=0.0;
        *yR=30.0;
        *NelemX=150;
        *NelemY=100;
        *PeriodicBD_X=0;
        *PeriodicBD_Y=0;


        break;

    }

    default:
    {
        *fixedDomain = 0 ;
        *fixedDisc = 0 ;
        break;
    }
    }


}



